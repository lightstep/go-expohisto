// Copyright The OpenTelemetry Authors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

package conversion // import "github.com/lightstep/go-expohisto/conversion"

import (
	"fmt"
	"math"
	"sort"

	"github.com/lightstep/go-expohisto/structure"
)

// genericHistogram represents any sort of histogram data structure.
type genericHistogram interface {
	// binStart is the lower boundary
	// of bucket with bucket ID index.
	bucketStart(index int) (float64, error)

	// binLimit is the upper boundary
	// of bucket with bucket ID index.
	bucketLimit(index int) (float64, error)

	// countInBin is the sample count
	// in bucket with bucket ID index.
	countInBucket(index int) (uint64, error)

	// numBuckets returns how many buckets (n)
	// are stored in this histogram. Those
	// buckets must have IDs [0, n).
	numBuckets() int
}

// bin is a part of a cumulative distribution.
type bin struct {
	index int32
	count float64
}

type Histogram = structure.Histogram[float64]

func FromGenericHistogram(
	input genericHistogram,
) (*Histogram, error) {

	// Pick a scale that gives good coverage of the output histogram
	// subject to limits in place.  Say we can support 1000 buckets
	// and the range is 0.5 (smallest) and 16 (largest).  Calculate
	// number of binary orders:
	//
	//   Log2(Largest / Smallest) = 5
	//
	// Scale 7 implies 128 buckets per order of magnitude, which
	// yields 640 buckets.  This logic has to be applied to both
	// the negative and positive range, if present.  (E.g., if the
	// negative range is -0.5 to -64, then the Log2() = 7, and
	// 7*128 < 1000 so both ranges fit 1000 buckets at scale 7.)
	var scale := pickOutputScale(negativeMin, negativeMax, positiveMin, positiveMax)

	// the output is constructed here.  Note the code in
	// ../structure isn't 100% ready for this.  that structure
	// expects to be created in an empty state with the maximum
	// scale (20) and automatically lowers scale as new
	// measurements arrive.  This structure, since we know the
	// output scale up front, could be just an slice of counts.
	// However, it's a slice of counts that has to be correctly
	// offset, so for now the pseudo-code calls this a *Histogram.
	var output *Histogram

	processBin := func(i int) error {
		bucketStart, err := input.bucketStart(i)
		if err != nil {
			return err
		}
		bucketLimit, err := input.bucketLimit(i)
		if err != nil {
			return err
		}
		if !fmath.IsFinite(bucketStart) || !fmath.IsFinite(bucketLimit) {
			return fmt.Errorf("bucket not finite: [%v, %v)", bucketStart, bucketLimit)
		}

		bucketWidth := bucketLimit - bucketStart
		if bucketStart > bucketLimit {
			return fmt.Errorf(
				"bucket has negative width: [%v, %v) [width=%v]",
				bucketStart,
				bucketLimit,
				bucketWidth,
			)
		}

		inputCount, err := input.countInBucket(i)
		if err != nil || count == 0 {
			return err
		}

		// entitledBin is a bin that knows the fractional
		// count it is entitled to.
		type entitledBin struct {
			index int32
			fractional float64
		}

		// remVals is the set of output buckets that have
		// fractional-count entitlement.
		var remVals []entitledBin

		// recordBin incre
		recordBin := func(
			idx int32,
			newCount uint64,
			fractional float64,
		) {
			remVals = append(remVals, entitledBin{
				index:      idx,
				fractional: fractional,
			})

			output.incrementIndex(idx, newCount)
		}

		if bucketStart <= 0 && bucketLimit >= 0 {
			// An input bucket that crosses zero spans a huge
			// number of output buckets.  Note: We
			// could split this into some number of bins
			// on either side of zero, within reason, but
			// this rule to create one bin is easiest.
			// Creating arbitrarily-small-exponent output
			// buckets near zero does not feel right.
			// Users should not do this (spam zero).
			recordBin(zeroBucket().binParams, inputCount, 0)
		} else {
			current := bucketStart
			currentIdx, err := lookupOutputBucketIndex(current)
			if err != nil {
				return fmt.Errorf("%w: [%v, %v)", err, bucketStart, bucketLimit)
			}

			for current < bucketLimit {
				outputStart := output.leftBoundary(currentIdx)
				outputLimit := output.rightBoundary(currentIdx)

				// Find what fraction of the original input interval
				// our new interval overlaps with.
				// Linearly interpolate into circllhist buckets.
				overlapWidth := (math.Min(circLimit, bucketLimit)) - (math.Max(circStart, bucketStart))
				density := float64(count) / bucketWidth

				fractionalCount := overlapWidth * density
				newCount := uint64(fractionalCount)

				recordBin(
					currentB,
					newCount,
					fractionalCount-float64(newCount),
				)

				// Done with this circllhist bucket.
				currentB, err = currentB.next()
				if err == errNextOverflow {
					// Special case recognizing at
					// some point there are no
					// more buckets.
					break
				}
				if err != nil {
					// Impossible, here.  OTLP
					// histograms cover the full
					// range.  If there are
					// output-range restrictions,
					// handle this error.
				}
				current = currentB.leftBoundary()
			}
		}

		if len(remVals) == 0 {
			// Impossible!
		}

		// Note: We expect that remainingCount is smaller than
		// len(remVals), i.e. that there was at most one unit of error
		// accumulated per bucket.
		remainingCount := inputCount - output.Count()

		// Sort by decreasing distance from how much
		// count each bucket is entitled to.
		sort.Slice(remVals, func(i, j int) bool {
			ei := remVals[i]
			ej := remVals[j]
			l := ei.fractional
			r := ej.fractional
			return l > r // Largest difference first.
		})

		// Apply the remaining count towards the buckets that
		// are missing the most of their entitled count.
		//
		// Every bucket gets perBucketBaseline. The first N buckets
		// get an additional one unit each.

		perBucketBaseline := remainingCount / uint64(len(remVals))
		remainingCount -= perBucketBaseline * uint64(len(remVals))

		// Pseudocode:
		// add perBucketBaseline to every item in `remVals`
		// add 1 to remainingCount remaining items in some order
		//
		// for ... {
		//    output.incrementIndex(..., perBucketBaseline)
		// }
		// for ... {
		//    output.incrementIndex(..., 1)
		// }
		//
		// Note there could definitely be a probabilistic
		// approach here or anywhere above!

		return nil
	}

	for i := 0; i < input.numBuckets(); i++ {
		if err := processBin(i); err != nil {
			return nil, err
		}
	}

	// Note: in an earlier version of this code there was another
	// stage of calculation here to format the output histogram
	// from an intermediate state.  Here, the pseudocode has the
	// output constructed on the fly, in pseudo-code, via an
	// unimplemented incrementIndex().  The code in ../structure
	// doesn't support that API and we don't need its full dynamic
	// structure to construct a histogram here.
	return output, nil
}

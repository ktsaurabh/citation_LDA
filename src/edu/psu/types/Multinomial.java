package edu.psu.types;



/**
* A probability distribution over a set of features represented as a {@link cc.mallet.types.FeatureVector}.
* The values associated with each element in the Multinomial/FeaturVector are probabilities
* and should sum to 1.
* Features are indexed using feature indices - the index into the underlying Alphabet -
* rather than using locations the way FeatureVectors do.
* <p>
* {@link cc.mallet.types.Multinomial.Estimator} provides a subhierachy
* of ways to generate an estimate of the probability distribution from counts associated
* with the features.
*
* @author Andrew McCallum <a href="mailto:mccallum@cs.umass.edu">mccallum@cs.umass.edu</a>
*/
public class Multinomial {
}

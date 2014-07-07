package edu.psu.topic;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.Vector;
import java.util.regex.Pattern;
import java.io.*;
import gnu.trove.*;

import edu.psu.types.*;
import edu.psu.util.*;

import com.panayotis.gnuplot.GNUPlotParameters;
import com.panayotis.gnuplot.JavaPlot;
import com.panayotis.gnuplot.dataset.FileDataSet;
import com.panayotis.gnuplot.layout.StripeLayout;
import com.panayotis.iodebug.Debug;
import com.panayotis.gnuplot.plot.AbstractPlot;
import com.panayotis.gnuplot.plot.DataSetPlot;
import com.panayotis.gnuplot.style.NamedPlotColor;
import com.panayotis.gnuplot.style.PlotStyle;
import com.panayotis.gnuplot.style.Style;
import com.panayotis.gnuplot.swing.JPlot;
import com.panayotis.gnuplot.terminal.PostscriptTerminal;
import com.panayotis.gnuplot.terminal.SVGTerminal;

public class linkPlsaLDA extends Model {

	/**
	 * @param args
	 */
	public int numSamples = 0;
	public int[] citedCountsPerTopic;
	public int citedCounts;
	public HashMap<Integer, Vector<Integer>> citedAssignments;

	private static LabelAlphabet newLabelAlphabet(int numTopics) {
		LabelAlphabet ret = new LabelAlphabet();
		for (int i = 0; i < numTopics; i++)
			ret.lookupIndex(i);
		return ret;
	}

	public Corpus readData(String directory) {
		String[] files = ReadDirectory.list(directory);
		System.out.println("There are total " + files.length + " Files in "
				+ directory);
		String line = null;
		Corpus corpus = new Corpus(files.length);

		for (int i = 0; i < files.length; i++) {
			String name = files[i].split("/")[files[i].split("/").length - 1]
					.replace(".txt", "");// check!!!
			if (!corpus.docAlphabet.contains(name)) {
				corpus.docAlphabet.lookupIndex(name);
			}
		}

		try {
			for (int i = 0; i < files.length; i++) {
				String name = files[i].split("/")[files[i].split("/").length - 1]
						.replace(".txt", "");// check!!!
				BufferedReader br = new BufferedReader(new FileReader(files[i]));
				try {

					Vector<String> lines = new Vector<String>();
					while ((line = br.readLine()) != null) {

						lines.add(line);
					}
					ContextDocument doc = new ContextDocument();
					doc.add(name,corpus.docAlphabet.lookupIndex(name), lines);
					corpus.addDocument(doc,
							corpus.docAlphabet.lookupIndex(name));
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}

			}
			System.out.println("Total Cited Documents:"
					+ corpus.citationAlphabet.size());
			System.out.println("Total Documents:" + corpus.docs.size());
			System.out.println("Total Vcabulary Size:"
					+ corpus.vocabulary.size());
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return corpus;

	}

	public void setWindowLength(Corpus corpus, int length) {
		for (int i = 0; i < corpus.docs.size(); i++) {
			ContextDocument doc = (ContextDocument) corpus.docs.get(i);
			for (int j = 1; j <= length; j++) {
				TIntObjectHashMap<Vector<Integer>> contextObject = doc.contextObject;
				int[] keys = contextObject.keys();
				for (int k = 0; k < keys.length; k++) {
					Vector<Integer> tmpContext = contextObject.get(keys[k]);
					if (keys[k] + j < doc.sentenseLength) {
						contextObject.put(keys[k] + j,
								(Vector<Integer>) tmpContext.clone());
					}
					if (keys[k] - j >= 0) {
						contextObject.put(keys[k] - j,
								(Vector<Integer>) tmpContext.clone());
					}
				}
			}
		}
	}

	public void InitializeParameters(Corpus corpus) {
		corpus.numTypes = corpus.vocabulary.size();

		corpus.beta = 0.001;
		corpus.betaSum = corpus.beta * corpus.numTypes;
		corpus.numCitations = corpus.docAlphabet.size();
		corpus.alpha = new double[numTopics];
		corpus.alphaSum = 0;
		for (int i = 0; i < numTopics; i++) {
			corpus.alpha[i] = 50.0 / (double) numTopics;
			corpus.alphaSum += corpus.alpha[i];
		}
		corpus.gamma = 0.001;
		corpus.gammaSum = corpus.gamma * corpus.numCitations;
		corpus.typeTopicCounts = new TIntIntHashMap[corpus.numTypes];
		corpus.tokensPerTopic = new int[numTopics];
		corpus.citationsPerTopic = new int[numTopics];
		corpus.phi_train = new double[corpus.vocabulary.size()][numTopics];
		corpus.theta_train = new double[corpus.docAlphabet.size()][numTopics];
		corpus.psi_train = new double[corpus.docAlphabet.size()][numTopics];
		corpus.citationTopicCounts = new TIntIntHashMap[corpus.numCitations];
		for (int i = 0; i < corpus.citationTopicCounts.length; i++)
			corpus.citationTopicCounts[i] = new TIntIntHashMap();
		for (int i = 0; i < corpus.typeTopicCounts.length; i++)
			corpus.typeTopicCounts[i] = new TIntIntHashMap();

	}

	public void InitializeCitedAssignments(Corpus corpus) {
		// initilize word, and citation factors
		this.random = new Randoms();
		citedAssignments = new HashMap<Integer, Vector<Integer>>();
		citedCountsPerTopic = new int[numTopics];
		// TODO AKM July 18: Why wasn't the next line there previously?
		// this.typeTopicCounts = newTypeTopicCounts;
		for (int doc_index = 0; doc_index < corpus.docs.size(); doc_index++) {
			if (Corpus.citationAlphabet.contains(doc_index)) {
				int[] words = corpus.docs.get(doc_index).wordSequence
						.getFeatures();
				Vector<Integer> topics = new Vector<Integer>();
				for (int j = 0; j < words.length; j++) {
					int topic = r.nextInt(numTopics);
					topics.add(topic);
					citedCountsPerTopic[topic]++;
					corpus.typeTopicCounts[words[j]].adjustOrPutValue(topic, 1,
							1);
					corpus.tokensPerTopic[topic]++;
					corpus.citationTopicCounts[doc_index].adjustOrPutValue(
							topic, 1, 1);
					corpus.citationsPerTopic[topic]++;
				}
				citedAssignments.put(doc_index, topics);
			}
		}

	}

	public void InitializeAssignments(Corpus corpus, LabelAlphabet topicAlphabet) {
		// initilize word, and citation factors
		this.random = new Randoms();

		// TODO AKM July 18: Why wasn't the next line there previously?
		// this.typeTopicCounts = newTypeTopicCounts;

		for (int doc_index = 0; doc_index < corpus.docs.size(); doc_index++) {
			ContextDocument doc = (ContextDocument) corpus.docs.get(doc_index);
			doc.topicAssignments = new LabelSequence(topicAlphabet,
					doc.wordSequence.getLength());

			int[] topics = doc.topicAssignments.getFeatures();
			int[] words = doc.wordSequence.getFeatures();
			// put a check to see if wordlength exceeds words.length
			for (int i = 0; i < doc.sentenseLength; i++) {
				for (int j = doc.SentenseBoundaries.get(i).firstElement(); j < doc.SentenseBoundaries
						.get(i).lastElement(); j++) {

					topics[j] = r.nextInt(numTopics);
					corpus.typeTopicCounts[words[j]].adjustOrPutValue(
							topics[j], 1, 1);
					corpus.tokensPerTopic[topics[j]]++;

				}
			}
			doc.citationTopicAssignment = new TIntIntHashMap();
			int[] keys = doc.citationSet.keys();
			for (int t = 0; t < keys.length; t++) {
				int topic = r.nextInt(numTopics);
				corpus.citationAlphabet.lookupIndex(keys[t], true);// initialize
																	// the
																	// citation
																	// alphabet
				doc.citationTopicAssignment.put(keys[t], topic);// put docid and
																// topic index
																// in
																// citationtopicassignment
				corpus.citationTopicCounts[keys[t]].adjustOrPutValue(topic, 1,
						1);
				corpus.citationsPerTopic[topic]++;
			}
		}

	}

	public void sampleCitedDocument(Corpus corpus, int docid,
			Vector<Integer> topicAssignments, int[] wordAssignments) {

		// decrement current sampling word

		// calculate document factor

		for (int i = 0; i < topicAssignments.size(); i++) {

			corpus.typeTopicCounts[wordAssignments[i]].adjustOrPutValue(
					topicAssignments.get(i), -1, 0);
			corpus.tokensPerTopic[topicAssignments.get(i)]--;
			corpus.citationTopicCounts[docid].adjustOrPutValue(
					topicAssignments.get(i), -1, 0);
			corpus.citationsPerTopic[topicAssignments.get(i)]--;
			citedCountsPerTopic[topicAssignments.get(i)]--;
			// doc_topics.adjustOrPutValue(topics[j], -1, 0);
			TIntIntHashMap currentTypeTopicCounts = corpus.typeTopicCounts[wordAssignments[i]];
			double[] topicDistribution = new double[numTopics];
			double topicDistributionSum = 0, weight = 0;
			for (int t = 0; t < numTopics; t++) {
				weight = ((currentTypeTopicCounts.get(t) + corpus.beta) / (corpus.tokensPerTopic[t] + corpus.betaSum))
						* ((corpus.citationTopicCounts[docid].get(t) + corpus.gamma) / (corpus.citationsPerTopic[t] + corpus.gammaSum))
						* (citedCountsPerTopic[t]);
				topicDistributionSum += weight;
				topicDistribution[t] = weight;
			}
			// System.out.print("Sampled:");
			// System.out.println(random.nextDiscrete(topicDistribution,
			// topicDistributionSum));
			int topic = random.nextDiscrete(topicDistribution,
					topicDistributionSum);
			topicAssignments.set(i, topic);
			// topics[j] = r.nextInt(numTopics);
			corpus.typeTopicCounts[wordAssignments[i]].adjustOrPutValue(topic,
					1, 1);

			corpus.tokensPerTopic[topic]++;
			corpus.citationTopicCounts[docid].adjustOrPutValue(topic, 1, 1);
			corpus.citationsPerTopic[topic]++;
			citedCountsPerTopic[topic]++;

		}

	}

	public void sampleOneDocument(Corpus corpus, ContextDocument doc) {
		TIntIntHashMap doc_topics = new TIntIntHashMap();
		// sample for each position
		int num_sentenses = doc.sentenseLength;
		int[] topics = doc.topicAssignments.getFeatures();
		int[] words = doc.wordSequence.getFeatures();
		for (int i = 0; i < topics.length; i++) {
			doc_topics.adjustOrPutValue(topics[i], 1, 1);
		}
		doc.citationTopicAssignment = new TIntIntHashMap();
		int[] keys = doc.citationTopicAssignment.keys();
		for (int t = 0; t < keys.length; t++) {
			doc_topics.adjustOrPutValue(
					doc.citationTopicAssignment.get(keys[t]), 1, 1);
		}

		for (int i = 0; i < doc.sentenseLength; i++) {
			for (int j = doc.SentenseBoundaries.get(i).firstElement(); j < doc.SentenseBoundaries
					.get(i).lastElement(); j++) {

				doc_topics.adjustOrPutValue(topics[j], -1, 0);
				corpus.typeTopicCounts[words[j]].adjustOrPutValue(topics[j],
						-1, 0);
				corpus.tokensPerTopic[topics[j]]--;
				// doc_topics.adjustOrPutValue(topics[j], -1, 0);
				TIntIntHashMap currentTypeTopicCounts = corpus.typeTopicCounts[words[j]];
				double[] topicDistribution = new double[numTopics];
				double topicDistributionSum = 0, weight = 0;
				for (int t = 0; t < numTopics; t++) {
					weight = ((currentTypeTopicCounts.get(t) + corpus.beta) / (corpus.tokensPerTopic[t] + corpus.betaSum))
							* ((doc_topics.get(t) + corpus.alpha[t]));
					topicDistributionSum += weight;
					topicDistribution[t] = weight;
				}
				// System.out.print("Sampled:");
				// System.out.println(random.nextDiscrete(topicDistribution,
				// topicDistributionSum));
				topics[j] = random.nextDiscrete(topicDistribution,
						topicDistributionSum);
				// topics[j] = r.nextInt(numTopics);
				corpus.typeTopicCounts[words[j]].adjustOrPutValue(topics[j], 1,
						1);
				doc_topics.adjustOrPutValue(topics[j], 1, 1);
				corpus.tokensPerTopic[topics[j]]++;
			}
		}
		keys = doc.citationTopicAssignment.keys();
		for (int c = 0; c < keys.length; c++) {
			int topic = doc.citationTopicAssignment.get(keys[c]);
			doc_topics.adjustOrPutValue(topic, -1, 0);
			corpus.citationTopicCounts[keys[c]].adjustOrPutValue(topic, -1, 0);
			corpus.citationsPerTopic[topic]--;
			TIntIntHashMap currentTypeTopicCounts = corpus.citationTopicCounts[keys[c]];
			double[] topicDistribution = new double[numTopics];
			double topicDistributionSum = 0, weight = 0;
			for (int t = 0; t < numTopics; t++) {
				weight = ((currentTypeTopicCounts.get(t) + corpus.gamma) / (corpus.citationsPerTopic[t] + corpus.gammaSum))
						* ((doc_topics.get(t) + corpus.alpha[t]));
				topicDistributionSum += weight;
				topicDistribution[t] = weight;
			}
			// System.out.print("Sampled:");
			// System.out.println(random.nextDiscrete(topicDistribution,
			// topicDistributionSum));
			topic = random
					.nextDiscrete(topicDistribution, topicDistributionSum);
			// topics[j] = r.nextInt(numTopics);
			corpus.citationTopicCounts[keys[c]].adjustOrPutValue(topic, 1, 1);
			doc_topics.adjustOrPutValue(topic, 1, 1);
			corpus.citationsPerTopic[topic]++;
		}

	}

	public void sampleCorpus(Corpus corpus, int iterations,
			boolean isContextAware) {
		
		this.InitializeCitedAssignments(corpus);
		
			for (int i = 0; i < corpus.docs.size(); i++) {
				sampleOneDocument(corpus, (ContextDocument) corpus.docs.get(i));
			}
			if (iterations % 10 == 0 && iterations>30) {
				for (int iter = 0; iter < 20; iter++) {
				for (int i = 0; i < corpus.docs.size(); i++) {
					if (Corpus.citationAlphabet.contains(i)) {
						this.sampleCitedDocument(corpus, i,
								citedAssignments.get(i),
								corpus.docs.get(i).wordSequence.getFeatures());
					}
				}
			}
			System.out.println("Iter=" + iterations);
		}
	}

	public void estimateParameters(Corpus corpus) {
		// estimate parameters after every k iterations after burn-in
		numSamples++;
		// int[][] typeTopicCounts = new
		// int[corpus.vocabulary.size()][numTopics];
		// int[][] docTopicCounts = new int[corpus.docs.size()][numTopics];
		// int[][] citationTopicCounts = new int[corpus.docs.size()][numTopics];

		for (int doc = 0; doc < corpus.docs.size(); doc++) {
			for (int i = 0; i < numTopics; i++) {
				if (numSamples > 1)
					corpus.theta_train[doc][i] *= (numSamples - 1);
				corpus.theta_train[doc][i] += (corpus.alpha[i] / (corpus.docs
						.get(doc).docLength + corpus.alphaSum));
			}
			int[] topics = corpus.docs.get(doc).topicAssignments.getFeatures();
			// System.out.println("=========");
			// System.out.println(topics.length);
			// System.out.println(corpus.docs.get(doc).topicAssignments.getLength());
			for (int i = 0; i < topics.length; i++) {
				corpus.theta_train[doc][topics[i]] += (1.0 / (corpus.docs
						.get(doc).docLength + corpus.alphaSum));
			}
			if (numSamples > 1) {
				for (int i = 0; i < numTopics; i++) {
					corpus.theta_train[doc][i] /= numSamples;
				}
			}
		}
		for (int k = 0; k < corpus.vocabulary.size(); k++) {
			for (int i = 0; i < numTopics; i++) {
				if (numSamples > 1)
					corpus.phi_train[k][i] *= (numSamples - 1);
				corpus.phi_train[k][i] += ((corpus.typeTopicCounts[k].get(i) + corpus.beta) / (corpus.tokensPerTopic[i] + corpus.betaSum));
				if (numSamples > 1)
					corpus.phi_train[k][i] /= numSamples;
			}
		}
		for (int k = 0; k < corpus.docAlphabet.size(); k++) {
			for (int i = 0; i < numTopics; i++) {
				if (numSamples > 1)
					corpus.psi_train[k][i] *= (numSamples - 1);
				corpus.psi_train[k][i] += ((corpus.citationTopicCounts[k]
						.get(i) + corpus.gamma) / (corpus.citationsPerTopic[i] + corpus.gammaSum));
				if (numSamples > 1)
					corpus.psi_train[k][i] /= (numSamples);
			}
		}

	}

	public double modelLogLikelihood(Corpus corpus) {
		double logLikelihood = 0.0;
		int nonZeroTopics;

		TIntIntHashMap[] typeTopicCounts = new TIntIntHashMap[corpus.numTypes];
		int[] tokensPerTopic = new int[numTopics];
		int[] citationsPerTopic = new int[numTopics];

		TIntIntHashMap[] citationTopicCounts = new TIntIntHashMap[corpus.numCitations];
		for (int i = 0; i < corpus.citationTopicCounts.length; i++)
			citationTopicCounts[i] = new TIntIntHashMap();
		for (int i = 0; i < corpus.typeTopicCounts.length; i++)
			typeTopicCounts[i] = new TIntIntHashMap();

		for (int doc_index = 0; doc_index < corpus.docs.size(); doc_index++) {
			ContextDocument doc = (ContextDocument) corpus.docs.get(doc_index);

			int[] topics = doc.topicAssignments.getFeatures();
			int[] words = doc.wordSequence.getFeatures();
			// put a check to see if wordlength exceeds words.length
			for (int i = 0; i < doc.sentenseLength; i++) {
				for (int j = doc.SentenseBoundaries.get(i).firstElement(); j < doc.SentenseBoundaries
						.get(i).lastElement(); j++) {

					if (doc.contextObject.contains(i)) {

						typeTopicCounts[words[j]].adjustOrPutValue(topics[j],
								1, 1);
						tokensPerTopic[topics[j]]++;
						// for (int citations = 0; citations < doc.contextObject
						// .get(i).size(); citations++) {

						// alphabet
						// topics[j] = r.nextInt(numTopics);
						// corpus.typeTopicCounts[words[j]].adjustOrPutValue(
						// topics[j], 1, 1);
						citationTopicCounts[doc.contextObject.get(i).get(0)]
								.adjustOrPutValue(topics[j], 1, 1);
						// corpus.tokensPerTopic[topics[j]]++;
						citationsPerTopic[topics[j]]++;
						// }
					} else {

						typeTopicCounts[words[j]].adjustOrPutValue(topics[j],
								1, 1);
						tokensPerTopic[topics[j]]++;
					}
				}
			}
		}

		// The likelihood of the model is a combination of a
		// Dirichlet-multinomial for the words in each topic
		// and a Dirichlet-multinomial for the topics in each
		// document.

		// The likelihood function of a dirichlet multinomial is
		// Gamma( sum_i alpha_i ) prod_i Gamma( alpha_i + N_i )
		// prod_i Gamma( alpha_i ) Gamma( sum_i (alpha_i + N_i) )

		// So the log likelihood is
		// logGamma ( sum_i alpha_i ) - logGamma ( sum_i (alpha_i + N_i) ) +
		// sum_i [ logGamma( alpha_i + N_i) - logGamma( alpha_i ) ]

		// Do the documents first

		int[] topicCounts = new int[numTopics];
		double[] topicLogGammas = new double[numTopics];
		int[] docTopics;
		int sampleSize = 0;

		for (int topic = 0; topic < numTopics; topic++) {
			topicLogGammas[topic] = Dirichlet
					.logGammaStirling(corpus.alpha[topic]);
		}

		for (int doc = 0; doc < corpus.docs.size(); doc++) {
			LabelSequence topicSequence = (LabelSequence) corpus.docs.get(doc).topicAssignments;

			docTopics = topicSequence.getFeatures();
			sampleSize += docTopics.length;

			for (int token = 0; token < docTopics.length; token++) {
				topicCounts[docTopics[token]]++;
			}

			for (int topic = 0; topic < numTopics; topic++) {
				if (topicCounts[topic] > 0) {

					logLikelihood += (Dirichlet
							.logGammaStirling(corpus.alpha[topic]
									+ topicCounts[topic]) - topicLogGammas[topic]);
				}
			}

			// subtract the (count + parameter) sum term
			logLikelihood -= Dirichlet.logGammaStirling(corpus.alphaSum
					+ docTopics.length);

			Arrays.fill(topicCounts, 0);
		}

		// add the parameter sum term
		logLikelihood += corpus.docs.size()
				* Dirichlet.logGammaStirling(corpus.alphaSum);

		// And the topics

		// Count the number of type-topic pairs
		int nonZeroTypeTopics = 0;

		for (int type = 0; type < corpus.vocabulary.size(); type++) {
			int[] usedTopics = typeTopicCounts[type].keys();

			for (int topic : usedTopics) {
				int count = typeTopicCounts[type].get(topic);
				if (count > 0) {
					nonZeroTypeTopics++;
					logLikelihood += Dirichlet.logGammaStirling(corpus.beta
							+ count);
				}
			}
		}

		for (int topic = 0; topic < numTopics; topic++) {

			logLikelihood -= Dirichlet
					.logGammaStirling((corpus.beta * numTopics)
							+ tokensPerTopic[topic]);
		}

		logLikelihood += (Dirichlet.logGammaStirling(corpus.beta * numTopics))
				- (Dirichlet.logGammaStirling(corpus.beta) * nonZeroTypeTopics);

		return Math.exp(-1 * logLikelihood / sampleSize);
	}

	public double testPerplexity(int numSamples, Corpus corpus) {
		double ll = 0;
		int sampleSize = 0;
		for (int doc = 0; doc < corpus.docs.size(); doc++) {
			ContextDocument document = (ContextDocument) corpus.docs.get(doc);
			sampleSize += document.wordLength;
			TIntIntHashMap word_counts = corpus.docs.get(doc).word_counts;
			int[] keys = word_counts.keys();
			double tmp_ll = 0;
			for (int i = 0; i < keys.length; i++) {
				for (int t = 0; t < numTopics; t++) {
					tmp_ll += corpus.theta_train[doc][t]
							* corpus.phi_train[keys[i]][t];
				}
				ll += (Math.log(tmp_ll) * word_counts.get(keys[i]));
			}

		}
		System.out.println("Samples=" + sampleSize);

		return Math.exp(-1 * ll / sampleSize);
	}

	public double sampleLikelihood(int numSamples, Corpus corpus) {
		double ll = 0;
		for (int doc = 0; doc < corpus.docs.size(); doc++) {

			TIntIntHashMap word_counts = corpus.docs.get(doc).word_counts;
			int[] keys = word_counts.keys();
			double tmp_ll = 0;
			for (int i = 0; i < keys.length; i++) {
				for (int t = 0; t < numTopics; t++) {
					tmp_ll += corpus.theta_train[doc][t]
							* corpus.phi_train[keys[i]][t];
				}
				ll += Math.log(tmp_ll * word_counts.get(keys[i]));
			}

		}
		return ll;
	}

	public double sampleLikelihood(int numSamples, Corpus corpus,
			TIntIntHashMap split) {
		double ll = 0;
		for (int doc = 0; doc < corpus.docs.size(); doc++) {
			if (!split.contains(doc)) {
				TIntIntHashMap word_counts = corpus.docs.get(doc).word_counts;
				int[] keys = word_counts.keys();
				double tmp_ll = 0;
				for (int i = 0; i < keys.length; i++) {
					for (int t = 0; t < numTopics; t++) {
						tmp_ll += corpus.theta_train[doc][t]
								* corpus.phi_train[keys[i]][t];
					}
					ll += Math.log(tmp_ll * word_counts.get(keys[i]));
				}
			}
		}
		return ll;
	}

	public double empiricalLikelihood(int numSamples, Corpus corpus) {
		double[][] likelihoods = new double[corpus.docs.size()][numSamples];
		double[] multinomial = new double[corpus.numTypes];
		double[] topicDistribution, currentSample, currentWeights;
		Dirichlet topicPrior = new Dirichlet(corpus.alpha);

		int sample, doc, topic, type, token, seqLen;
		FeatureSequence fs;

		for (sample = 0; sample < numSamples; sample++) {
			topicDistribution = topicPrior.nextDistribution();
			Arrays.fill(multinomial, 0.0);

			for (topic = 0; topic < numTopics; topic++) {
				for (type = 0; type < corpus.numTypes; type++) {
					multinomial[type] += topicDistribution[topic]
							* (corpus.beta + corpus.typeTopicCounts[type]
									.get(topic))
							/ (corpus.betaSum + corpus.tokensPerTopic[topic]);
				}
			}

			// Convert to log probabilities
			for (type = 0; type < corpus.numTypes; type++) {
				assert (multinomial[type] > 0.0);
				multinomial[type] = Math.log(multinomial[type]);
			}

			for (doc = 0; doc < corpus.docs.size(); doc++) {
				fs = (FeatureSequence) corpus.docs.get(doc).wordSequence;
				seqLen = fs.getLength();

				for (token = 0; token < seqLen; token++) {
					type = fs.getIndexAtPosition(token);

					// Adding this check since testing instances may
					// have types not found in training instances,
					// as pointed out by Steven Bethard.
					if (type < corpus.numTypes) {
						likelihoods[doc][sample] += multinomial[type];
					}
				}
			}
		}

		double averageLogLikelihood = 0.0;
		double logNumSamples = Math.log(numSamples);
		for (doc = 0; doc < corpus.docs.size(); doc++) {
			double max = Double.NEGATIVE_INFINITY;
			for (sample = 0; sample < numSamples; sample++) {
				if (likelihoods[doc][sample] > max) {
					max = likelihoods[doc][sample];
				}
			}

			double sum = 0.0;
			for (sample = 0; sample < numSamples; sample++) {
				sum += Math.exp(likelihoods[doc][sample] - max);
			}

			averageLogLikelihood += Math.log(sum) + max - logNumSamples;
		}

		return averageLogLikelihood;

	}

	public void precision(Corpus corpus) {

	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		// Read documents
		String input = args[0];
		String output = args[1];
		int windowLength = Integer.parseInt(args[2]);
		int iterations = Integer.parseInt(args[3]);
		double[][] plot_citing = new double[windowLength][2];
		double[][] plot_cited = new double[windowLength][2];
		double[][] plot_citing_cited = new double[windowLength][2];
		double[] citingavg = new double[windowLength];
		double[] citedavg = new double[windowLength];
		double[] citing_cited_avg = new double[windowLength];
		for (int i = 0; i < windowLength; i++) {
			citing_cited_avg[i] = 0;
		}
		linkPlsaLDA linkPlsaLDA = new linkPlsaLDA();
		System.out.println("Reading Data.....");
		Corpus corpus = linkPlsaLDA.readData(input);
		System.out.println("Done");
		double[][][] citingKL = new double[windowLength][corpus.docs.size()][2];
		double[][][] citedKL = new double[windowLength][corpus.docs.size()][2];
		double[][][] citingcitedKL = new double[windowLength][corpus.docs
				.size()][2];
		for (int length = 0; length < windowLength; length++) {
			for (int i = 0; i < corpus.docs.size(); i++) {
				citingKL[length][i][0] = i;
				citingKL[length][i][1] = 0;
				citedKL[length][i][0] = i;
				citedKL[length][i][1] = 0;
				citingcitedKL[length][i][0] = i;
				citingcitedKL[length][i][1] = 0;
			}

		}
		for (int length = 0; length < windowLength; length++) {

			linkPlsaLDA = new linkPlsaLDA();
			System.out.println("Reading Data.....");
			corpus = linkPlsaLDA.readData(input);
			System.out.println("Done");
			System.out.println("Setting Window Length.....");
			linkPlsaLDA.setWindowLength(corpus, length);
			System.out.println("Done");
			System.out.println("Size of docs=" + corpus.docs.size());
			linkPlsaLDA.numTopics = 10;
			linkPlsaLDA.sampleCorpus(corpus, iterations, true);
			if (false) {

				int count = 0;
				for (int i = 0; i < corpus.docs.size(); i++) {
					ContextDocument doc = (ContextDocument) corpus.docs.get(i);
					int[] citations = doc.citationSet.keys();
					for (int j = 0; j < citations.length; j++) {
						count++;
					}

				}

				count = 0;
				for (int i = 0; i < corpus.docs.size(); i++) {
					ContextDocument doc = (ContextDocument) corpus.docs.get(i);
					int[] citations = doc.citationSet.keys();
					if (citations.length > 0) {
						count++;

					}
				}
				System.out.println(count);
				for (int i = 0; i < corpus.docs.size(); i++) {
					ContextDocument doc = (ContextDocument) corpus.docs.get(i);
					int[] citations = doc.citationSet.keys();
					for (int j = 0; j < citations.length; j++) {
						double[] tmpKL = new double[3];// linkPlsaLDA.calculateJSDiv(corpus,
														// i, citations[j],
														// shortc);
						if (tmpKL[0] > 5 || tmpKL[1] > 5 || tmpKL[2] > 5) {
							tmpKL[0] = 0;
							tmpKL[1] = 0;
							tmpKL[2] = 0;
						}
						citing_cited_avg[length] += (tmpKL[0] / count);
						citingavg[length] += (tmpKL[1] / count);
						citedavg[length] += (tmpKL[2] / count);

						// System.out.println(i + ":" + citations[j] + ":" +
						// tmpKL[0]
						// + ":" + tmpKL[1] + ":" + tmpKL[2]);
						citingcitedKL[length][i][0] = i;
						citingcitedKL[length][i][1] = tmpKL[0];
						citingKL[length][i][0] = i;
						citingKL[length][i][1] += (tmpKL[1] / citations.length);
						citedKL[length][i][0] = i;
						citedKL[length][i][1] += (tmpKL[2] / citations.length);
					}
				}
			}

		}
		System.exit(0);
		for (int length = 0; length < windowLength; length++) {
			System.out.println("Length=" + length + " Citing Avg="
					+ citingavg[length] + " Cited Avg=" + citedavg[length]);
			plot_citing[length][0] = length;
			plot_citing[length][1] = citingavg[length];
			plot_cited[length][0] = length;
			plot_cited[length][1] = citedavg[length];
			plot_citing_cited[length][0] = length;
			plot_citing_cited[length][1] = citing_cited_avg[length];

		}

		JavaPlot p;
		PostscriptTerminal epsf;
		DataSetPlot s;
		StripeLayout lo;
		String GNUPLOT = "/usr/local/bin/gnuplot";
		// for (int i = 0; i < windowLength; i++) {
		p = new JavaPlot(GNUPLOT);
		epsf = new PostscriptTerminal(output + "/plot_sc.eps");
		epsf.setColor(true);
		p.setTerminal(epsf);

		p.setTitle("");
		p.getAxis("x").setLabel("Citation context radius", "Arial", 20);
		p.getAxis("y").setLabel("JS Divegence", "Arial", 20);

		p.setKey(JavaPlot.Key.TOP_RIGHT);

		s = new DataSetPlot(plot_citing);
		p.addPlot(s);
		s = new DataSetPlot(plot_cited);
		p.addPlot(s);
		// s = new DataSetPlot(plot_citing_cited);
		// p.addPlot(s);

		((AbstractPlot) p.getPlots().get(0))
				.setTitle("Mean Divergence between Citing and Context");
		((AbstractPlot) p.getPlots().get(1))
				.setTitle("Mean Divergence between Cited and Context");
		// ((AbstractPlot) p.getPlots().get(2)).setTitle("Citing_Cited");
		PlotStyle stl = ((AbstractPlot) p.getPlots().get(0)).getPlotStyle();
		stl.setStyle(Style.LINESPOINTS);
		stl = ((AbstractPlot) p.getPlots().get(1)).getPlotStyle();
		stl.setStyle(Style.LINESPOINTS);
		// stl = ((AbstractPlot) p.getPlots().get(2)).getPlotStyle();
		// stl.setStyle(Style.LINESPOINTS);

		lo = new StripeLayout();
		lo.setColumns(99);
		p.getPage().setLayout(lo);
		p.plot();
		// }

	}

}

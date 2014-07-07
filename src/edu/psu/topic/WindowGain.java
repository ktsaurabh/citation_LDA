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

public class WindowGain {

	/**
	 * @param args
	 */
	public int numTopics;
	public Random r = new Random();
	public Randoms random;
	public static int numSamples = 0;

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

	public double[] calculateKLDiv(Corpus corpus, int doc, int citation) {
		// Given a citing document and cited document objects, it calculates
		// the KL divergence between citing document and context of cited
		// document in KL[0]
		// and the KL divergence between cited document and context of cited
		// document in KL[1]
		double[] KL = new double[2];
		int contextLength = 0, citingLength = 0, citedLength = 0;
		TIntIntHashMap citationWordCount = new TIntIntHashMap();
		TIntIntHashMap citingWordCount = new TIntIntHashMap();
		ContextDocument citing = (ContextDocument) corpus.docs.get(doc);
		ContextDocument cited = (ContextDocument) corpus.docs.get(citation);
		// populate citation context frequency vector
		int[] words = citing.wordSequence.getFeatures();
		// System.out.println(words.length);
		for (int i = 0; i < citing.sentenseLength; i++) {
			if (citing.contextObject.contains(i)) {
				// System.out.println(citing.SentenseBoundaries.get(i).firstElement()+":"+citing.SentenseBoundaries.get(i).lastElement());
				for (int j = citing.SentenseBoundaries.get(i).firstElement(); j < citing.SentenseBoundaries
						.get(i).lastElement(); j++) {

					if (citationWordCount.contains(words[j])) {
						citationWordCount.put(words[j], 1);
					} else {
						citationWordCount.put(words[j],
								citationWordCount.get(words[j]) + 1);
					}
					contextLength++;
				}

			}

		}
		// now calculate KL between citing and context
		double divergence = 0;
		words = citing.word_counts.keys();
		for (int i = 0; i < words.length; i++) {
			if (citing.word_counts.get(words[i])
					- citationWordCount.get(words[i]) > 0)
				citingWordCount.put(words[i], citing.word_counts.get(words[i])
						- citationWordCount.get(words[i]));
		}
		citingLength = citing.wordLength - contextLength;
		if (citingLength == 0) {
			KL[0] = 0;
			KL[1] = 0;
			return KL;
		}
		words = citingWordCount.keys();
		for (int i = 0; i < words.length; i++) {
			divergence += citingWordCount.get(words[i])
					* Math.log(citingWordCount.get(words[i]));
			if (citationWordCount.contains(words[i])) {
				divergence -= citingWordCount.get(words[i])
						* Math.log(citationWordCount.get(words[i]));
			}
		}
		System.out.println("citingLength=" + citingLength);
		System.out.println("contextLength=" + contextLength);
		divergence /= citingLength;
		divergence += Math.log(contextLength);
		divergence -= Math.log(citingLength);
		KL[0] = divergence;

		// now calculate KL between cited and context
		divergence = 0;
		words = cited.word_counts.keys();
		for (int i = 0; i < words.length; i++) {
			divergence += cited.word_counts.get(words[i])
					* Math.log(cited.word_counts.get(words[i]));
			if (citationWordCount.contains(words[i])) {
				divergence -= cited.word_counts.get(words[i])
						* Math.log(citationWordCount.get(words[i]));
			}
		}
		divergence /= cited.wordLength;
		divergence += Math.log(contextLength);
		divergence -= Math.log(cited.wordLength);
		KL[1] = divergence;
		return KL;
	}

	public double[] calculateJSDiv_0(Corpus corpus, int doc, int citation) {
		// Given a citing document and cited document objects, it calculates
		// the KL divergence between citing document and context of cited
		// document in KL[0]
		// and the KL divergence between cited document and context of cited
		// document in KL[1]
		double[] KL = new double[3];
		int contextLength = 0, citingLength = 0, citedLength = 0;
		TIntIntHashMap citationWordCount = new TIntIntHashMap();
		TIntIntHashMap citingWordCount = new TIntIntHashMap();
		ContextDocument citing = (ContextDocument) corpus.docs.get(doc);
		ContextDocument cited = (ContextDocument) corpus.docs.get(citation);
		// calculate citing and cited information difference
		citingWordCount = citing.word_counts;
		int[] words = citingWordCount.keys();
		citationWordCount = cited.word_counts;

		citingLength = citing.wordLength;
		contextLength = cited.wordLength;
		if (citingLength == 0 || contextLength == 0) {
			KL[0] = 0;
			KL[1] = 0;
			KL[2] = 0;
			return KL;
		}
		double divergence = 0;
		for (int i = 0; i < words.length; i++) {
			if (citationWordCount.contains(words[i])) {
				divergence += (((double) citingWordCount.get(words[i]) / (double) citingLength) * Math
						.log(((double) 2 * citingWordCount.get(words[i]) * contextLength)
								/ (double) (citingWordCount.get(words[i])
										* contextLength + citationWordCount
										.get(words[i]) * citingLength)));
			}
		}
		words = citationWordCount.keys();

		for (int i = 0; i < words.length; i++) {
			if (citingWordCount.contains(words[i])) {
				divergence += (((double) citationWordCount.get(words[i]) / (double) contextLength) * Math
						.log(((double) 2 * citationWordCount.get(words[i]) * citingLength)
								/ (double) (citationWordCount.get(words[i])
										* citingLength + citingWordCount
										.get(words[i]) * contextLength)));
			}
		}

		KL[0] = divergence / 2;

		return KL;
	}

	public double[] calculateJSDiv(Corpus corpus, int doc, int citation, int scn) {
		// Given a citing document and cited document objects, it calculates
		// the KL divergence between citing document and context of cited
		// document in KL[0]
		// and the KL divergence between cited document and context of cited
		// document in KL[1]
		double[] KL = new double[3];
		int contextLength = 0, citingLength = 0, citedLength = 0;
		TIntIntHashMap citationWordCount = new TIntIntHashMap();
		TIntIntHashMap citingWordCount = new TIntIntHashMap();
		ContextDocument citing = (ContextDocument) corpus.docs.get(doc);
		ContextDocument cited = (ContextDocument) corpus.docs.get(citation);
		// calculate citing and cited information difference
		citingWordCount = citing.word_counts;
		int[] words = citingWordCount.keys();
		citationWordCount = cited.word_counts;

		citingLength = citing.wordLength;
		contextLength = cited.wordLength;
		if (citingLength == 0 || contextLength == 0) {
			KL[0] = 0;
			KL[1] = 0;
			KL[2] = 0;
			return KL;
		}
		double divergence = 0;
		for (int i = 0; i < words.length; i++) {
			if (citationWordCount.contains(words[i])) {
				divergence += (((double) citingWordCount.get(words[i]) / (double) citingLength) * Math
						.log(((double) 2 * citingWordCount.get(words[i]) * contextLength)
								/ (double) (citingWordCount.get(words[i])
										* contextLength + citationWordCount
										.get(words[i]) * citingLength)))
						/ Math.log(2.0);
				;
			}
		}
		words = citationWordCount.keys();

		for (int i = 0; i < words.length; i++) {
			if (citingWordCount.contains(words[i])) {
				divergence += (((double) citationWordCount.get(words[i]) / (double) contextLength) * Math
						.log(((double) 2 * citationWordCount.get(words[i]) * citingLength)
								/ (double) (citationWordCount.get(words[i])
										* citingLength + citingWordCount
										.get(words[i]) * contextLength)))
						/ Math.log(2.0);
				;
			}
		}

		KL[0] = divergence / 2;
		citingWordCount = new TIntIntHashMap();
		citationWordCount = new TIntIntHashMap();
		contextLength = 0;

		// populate citation context frequency vector
		words = citing.wordSequence.getFeatures();
		// System.out.println(words.length);
		for (int i = 0; i < citing.sentenseLength; i++) {
			if (citing.contextObject.contains(i)) {
				// System.out.println(citing.SentenseBoundaries.get(i).firstElement()+":"+citing.SentenseBoundaries.get(i).lastElement());
				for (int j = citing.SentenseBoundaries.get(i).firstElement(); j < citing.SentenseBoundaries
						.get(i).lastElement(); j++) {

					if (!citationWordCount.contains(words[j])) {
						citationWordCount.put(words[j], 1);
					} else {
						citationWordCount.put(words[j],
								citationWordCount.get(words[j]) + 1);
					}
					contextLength++;
				}

			}

		}
		// now calculate KL between citing and context
		divergence = 0;
		words = citing.word_counts.keys();
		for (int i = 0; i < words.length; i++) {
			if (citing.word_counts.get(words[i])
					- citationWordCount.get(words[i]) > 0)
				citingWordCount.put(words[i], citing.word_counts.get(words[i])
						- citationWordCount.get(words[i]));
		}
		citingLength = citing.wordLength - contextLength;
		/*
		 * System.out.println("citing:" + corpus.docs.get(doc).docName +
		 * " citation:" + citation + " citingLength:" + citingLength +
		 * " contextLength:" + contextLength);
		 */
		if (citingLength == 0) {
			KL[1] = 0;
			KL[2] = 0;
			return KL;
		}
		if (scn == 1) {
			citingWordCount = citing.word_counts;// short-circuit
			citingLength = citing.wordLength;
		}
		words = citingWordCount.keys();
		for (int i = 0; i < words.length; i++) {
			if (citationWordCount.contains(words[i])) {
				divergence += (((double) citingWordCount.get(words[i]) / (double) citingLength) * Math
						.log(((double) 2 * citingWordCount.get(words[i]) * contextLength)
								/ (double) (citingWordCount.get(words[i])
										* contextLength + citationWordCount
										.get(words[i]) * citingLength)))
						/ Math.log(2.0);
				;
			}
		}
		words = citationWordCount.keys();
		for (int i = 0; i < words.length; i++) {
			if (citingWordCount.contains(words[i])) {
				divergence += (((double) citationWordCount.get(words[i]) / (double) contextLength) * Math
						.log(((double) 2 * citationWordCount.get(words[i]) * citingLength)
								/ (double) (citationWordCount.get(words[i])
										* citingLength + citingWordCount
										.get(words[i]) * contextLength)))
						/ Math.log(2.0);
				;
			}
		}

		KL[1] = divergence / 2;

		// now calculate KL between cited and context
		divergence = 0;
		words = cited.word_counts.keys();
		citedLength = cited.wordLength;
		if (citedLength == 0) {
			KL[1] = 0;
			KL[2] = 0;
			return KL;
		}
		for (int i = 0; i < words.length; i++) {
			if (citationWordCount.contains(words[i])) {
				divergence += ((double) cited.word_counts.get(words[i]) / (double) citedLength)
						* Math.log(((double) 2
								* cited.word_counts.get(words[i]) * contextLength)
								/ (double) (cited.word_counts.get(words[i])
										* contextLength + citationWordCount
										.get(words[i]) * citedLength))
						/ Math.log(2.0);
			}
		}
		words = citationWordCount.keys();
		for (int i = 0; i < words.length; i++) {
			if (cited.word_counts.contains(words[i])) {
				divergence += ((double) citationWordCount.get(words[i]) / (double) contextLength)
						* Math.log(((double) 2
								* citationWordCount.get(words[i]) * citedLength)
								/ (double) (citationWordCount.get(words[i])
										* citedLength + cited.word_counts
										.get(words[i]) * contextLength))
						/ Math.log(2.0);
			}
		}

		KL[2] = divergence / 2;
		return KL;
	}

	public void Initialize(Corpus corpus, LabelAlphabet topicAlphabet) {
		// initilize word, and citation factors

		corpus.numTypes = corpus.vocabulary.size();
		corpus.typeTopicCounts = new TIntIntHashMap[corpus.numTypes];
		corpus.tokensPerTopic = new int[numTopics];
		corpus.citationsPerTopic = new int[numTopics];
		corpus.phi_train = new double[corpus.vocabulary.size()][numTopics];
		corpus.theta_train = new double[corpus.docs.size()][numTopics];
		corpus.psi_train = new double[corpus.citationAlphabet.size()][numTopics];
		for (int i = 0; i < corpus.typeTopicCounts.length; i++)
			corpus.typeTopicCounts[i] = new TIntIntHashMap();
		// TODO AKM July 18: Why wasn't the next line there previously?
		// this.typeTopicCounts = newTypeTopicCounts;
		corpus.betaSum = corpus.beta * corpus.numTypes;

		corpus.numCitations = corpus.citationAlphabet.size();
		corpus.citationTopicCounts = new TIntIntHashMap[corpus.numCitations];
		for (int i = 0; i < corpus.citationTopicCounts.length; i++)
			corpus.citationTopicCounts[i] = new TIntIntHashMap();
		// TODO AKM July 18: Why wasn't the next line there previously?
		// this.typeTopicCounts = newTypeTopicCounts;
		corpus.gammaSum = corpus.gamma * corpus.numCitations;

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
					if (doc.contextObject.contains(i)) {
						topics[j] = r.nextInt(numTopics);
						corpus.typeTopicCounts[words[j]].adjustOrPutValue(
								topics[j], 1, 1);
						corpus.citationTopicCounts[doc.contextObject.get(i)
								.firstElement()].adjustOrPutValue(topics[j], 1,
								1);
						corpus.tokensPerTopic[topics[j]]++;
						corpus.citationsPerTopic[topics[j]]++;
					} else {
						topics[j] = r.nextInt(numTopics);
						corpus.typeTopicCounts[words[j]].adjustOrPutValue(
								topics[j], 1, 1);
						corpus.tokensPerTopic[topics[j]]++;
					}
				}

			}

			/*
			 * for (int i = 0; i < doc.wordSequence.getLength(); i++) {
			 * topics[i] = r.nextInt();
			 * corpus.typeTopicCounts[words[i]].adjustOrPutValue(topics[i], 1,
			 * 1);
			 * 
			 * }
			 */

		}
	}

	public void sampleOneDocument(Corpus corpus, ContextDocument doc) {

		// decrement current sampling word

		// calculate document factor
		TIntIntHashMap doc_topics = new TIntIntHashMap();
		// sample for each position
		int num_sentenses = doc.sentenseLength;
		int[] topics = doc.topicAssignments.getFeatures();
		int[] words = doc.wordSequence.getFeatures();
		for (int i = 0; i < topics.length; i++) {
			doc_topics.adjustOrPutValue(topics[i], 1, 1);
		}

		for (int i = 0; i < doc.sentenseLength; i++) {
			for (int j = doc.SentenseBoundaries.get(i).firstElement(); j < doc.SentenseBoundaries
					.get(i).lastElement(); j++) {
				if (doc.contextObject.contains(i)) {
					doc_topics.adjustOrPutValue(topics[j], -1, 0);// check for
																	// the
																	// adjust
																	// and put
																	// order
					corpus.typeTopicCounts[words[j]].adjustOrPutValue(
							topics[j], -1, 0);
					corpus.citationTopicCounts[doc.contextObject.get(i)
							.firstElement()].adjustOrPutValue(topics[j], -1, 0);
					corpus.tokensPerTopic[topics[j]]--;
					corpus.citationsPerTopic[topics[j]]--;
					doc_topics.adjustOrPutValue(topics[j], -1, 0);
					TIntIntHashMap currentTypeTopicCounts = corpus.typeTopicCounts[words[j]];
					TIntIntHashMap currentCitationTopicCounts = corpus.citationTopicCounts[doc.contextObject
							.get(i).firstElement()];
					double[] topicDistribution = new double[numTopics];
					double topicDistributionSum = 0, weight = 0;
					for (int t = 0; t < numTopics; t++) {
						weight = ((currentTypeTopicCounts.get(t) + corpus.beta) / (corpus.tokensPerTopic[t] + corpus.betaSum))
								* ((currentCitationTopicCounts.get(t) + corpus.gamma) / (corpus.citationsPerTopic[t] + corpus.gammaSum))
								* ((doc_topics.get(t) + corpus.alpha[t]));
						topicDistributionSum += weight;
						topicDistribution[t] = weight;
					}
					topics[j] = random.nextDiscrete(topicDistribution,
							topicDistributionSum);
					// topics[j] = r.nextInt(numTopics);
					corpus.typeTopicCounts[words[j]].adjustOrPutValue(
							topics[j], 1, 1);
					corpus.citationTopicCounts[doc.contextObject.get(i)
							.firstElement()].adjustOrPutValue(topics[j], 1, 1);
					doc_topics.adjustOrPutValue(topics[j], 1, 1);
					corpus.tokensPerTopic[topics[j]]++;
					corpus.citationsPerTopic[topics[j]]++;
				} else {
					doc_topics.adjustOrPutValue(topics[j], -1, 0);// check for
																	// the
																	// adjust
																	// and put
																	// order
					corpus.typeTopicCounts[words[j]].adjustOrPutValue(
							topics[j], -1, 0);
					corpus.citationTopicCounts[doc.contextObject.get(i)
							.firstElement()].adjustOrPutValue(topics[j], -1, 0);
					corpus.tokensPerTopic[topics[j]]--;
					corpus.citationsPerTopic[topics[j]]--;
					doc_topics.adjustOrPutValue(topics[j], -1, 0);
					TIntIntHashMap currentTypeTopicCounts = corpus.typeTopicCounts[words[j]];
					double[] topicDistribution = new double[numTopics];
					double topicDistributionSum = 0, weight = 0;
					for (int t = 0; t < numTopics; t++) {
						weight = ((currentTypeTopicCounts.get(t) + corpus.beta) / (corpus.tokensPerTopic[t] + corpus.betaSum))
								* ((doc_topics.get(t) + corpus.alpha[t]));
						topicDistributionSum += weight;
						topicDistribution[t] = weight;
					}
					topics[j] = random.nextDiscrete(topicDistribution,
							topicDistributionSum);
					// topics[j] = r.nextInt(numTopics);
					corpus.typeTopicCounts[words[j]].adjustOrPutValue(
							topics[j], 1, 1);
					corpus.citationTopicCounts[doc.contextObject.get(i)
							.firstElement()].adjustOrPutValue(topics[j], 1, 1);
					doc_topics.adjustOrPutValue(topics[j], 1, 1);
					corpus.tokensPerTopic[topics[j]]++;
					corpus.citationsPerTopic[topics[j]]++;
				}
			}

		}

	}

	public void sampleCorpus(Corpus corpus) {
		for (int i = 0; i < corpus.docs.size(); i++) {
			sampleOneDocument(corpus, (ContextDocument) corpus.docs.get(i));
		}
	}

	public void estimateParameters(Corpus corpus) {
		numSamples++;
		int[][] typeTopicCounts = new int[corpus.vocabulary.size()][numTopics];
		int[][] docTopicCounts = new int[corpus.docs.size()][numTopics];
		int[][] citationTopicCounts = new int[corpus.docs.size()][numTopics];
		for (int doc = 0; doc < corpus.docs.size(); doc++) {
			for (int i = 0; i < numTopics; i++) {
				corpus.theta_train[doc][i] = corpus.alpha[i];
			}
			int[] topics = corpus.docs.get(doc).topicAssignments.getFeatures();
			for (int i = 0; i < topics.length; i++) {
				corpus.theta_train[doc][topics[i]] += 1;
			}
			for (int i = 0; i < numTopics; i++) {
				corpus.theta_train[doc][topics[i]] /= (corpus.docs.get(doc).docLength + corpus.alphaSum);
			}
		}

	}

	public double sampleLikelihood(int numSamples, Corpus corpus) {
		double ll = 0;
		for (int doc = 0; doc < corpus.docs.size(); doc++) {
			TIntIntHashMap word_counts = corpus.docs.get(doc).word_counts;

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
		int shortc = Integer.parseInt(args[3]);
		double[][] plot_citing = new double[windowLength][2];
		double[][] plot_cited = new double[windowLength][2];
		double[][] plot_citing_cited = new double[windowLength][2];
		double[] citingavg = new double[windowLength];
		double[] citedavg = new double[windowLength];
		double[] citing_cited_avg = new double[windowLength];
		for (int i = 0; i < windowLength; i++) {
			citing_cited_avg[i] = 0;
		}
		citeLDA citelda = new citeLDA();
		System.out.println("Reading Data.....");
		Corpus corpus = citelda.readData(input);
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

			citelda = new citeLDA();
			System.out.println("Reading Data.....");
			corpus = citelda.readData(input);
			System.out.println("Done");
			System.out.println("Setting Window Length.....");
			citelda.setWindowLength(corpus, length);
			System.out.println("Done");
			System.out.println("Size of docs=" + corpus.docs.size());

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
					double[] tmpKL = null;//citelda.calculateJSDiv(corpus, i,
							//citations[j], shortc);
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
		if (shortc == 0)// extract context out from citing and then do the comparision
			epsf = new PostscriptTerminal(output + "/plot_scn.eps");
		else
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
		//s = new DataSetPlot(plot_citing_cited);
		//p.addPlot(s);

		((AbstractPlot) p.getPlots().get(0)).setTitle("Mean Divergence between Citing and Context");
		((AbstractPlot) p.getPlots().get(1)).setTitle("Mean Divergence between Cited and Context");
		//((AbstractPlot) p.getPlots().get(2)).setTitle("Citing_Cited");
		PlotStyle stl = ((AbstractPlot) p.getPlots().get(0)).getPlotStyle();
		stl.setStyle(Style.LINESPOINTS);
		stl = ((AbstractPlot) p.getPlots().get(1)).getPlotStyle();
		stl.setStyle(Style.LINESPOINTS);
		//stl = ((AbstractPlot) p.getPlots().get(2)).getPlotStyle();
		//stl.setStyle(Style.LINESPOINTS);

		lo = new StripeLayout();
		lo.setColumns(99);
		p.getPage().setLayout(lo);
		p.plot();
		// }

	}

}

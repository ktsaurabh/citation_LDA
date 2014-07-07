package edu.psu.factorization;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Random;
import java.util.TreeSet;
import java.util.Vector;
import java.util.regex.Pattern;
import java.io.*;

import gnu.trove.*;

import edu.psu.topic.Model;
import edu.psu.topic.citeLDA;
import edu.psu.topic.linkLDA;
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

public class NMF {

	/**
	 * @param args
	 */
	public int numTopics;
	public static Random r = new Random();
	public Randoms random;
	public double[][] citationFactors;
	public double[][] wordFactors;
	public double[][] documentFactors;
	public double[][] ApproximatedValue;
	public double[][] ApproximatedWordValue;
	public double[][] interDocFactors;
	public double[][] interCitationFactors;
	public double[][] interWordFactors;
	public double[][] interDocWordFactors;
	public double[] docFactorSum;
	public double[] citationFactorSum;
	public double[] wordFactorSum;

	public void initializeFactors(Corpus corpus, int factors) {
		citationFactors = new double[corpus.docs.size()][factors];
		interCitationFactors = new double[corpus.docs.size()][factors];
		wordFactors = new double[corpus.vocabulary.size()][factors];
		interWordFactors = new double[corpus.vocabulary.size()][factors];
		documentFactors = new double[corpus.docs.size()][factors];
		interDocFactors = new double[corpus.docs.size()][factors];
		interDocWordFactors = new double[corpus.docs.size()][factors];
		ApproximatedValue = new double[corpus.docs.size()][corpus.docs.size()];
		ApproximatedWordValue = new double[corpus.docs.size()][corpus.vocabulary
				.size()];
		docFactorSum = new double[factors];
		citationFactorSum = new double[factors];
		wordFactorSum = new double[factors];

		for (int doc = 0; doc < corpus.docs.size(); doc++) {
			for (int factor = 0; factor < factors; factor++) {
				documentFactors[doc][factor] = r.nextDouble();
			}
		}
		for (int citation = 0; citation < corpus.docs.size(); citation++) {
			if (corpus.citationAlphabet.contains(citation)) {
				for (int factor = 0; factor < factors; factor++) {
					citationFactors[citation][factor] = r.nextDouble();
				}
			}
		}
		for (int word = 0; word < corpus.vocabulary.size(); word++) {
			for (int factor = 0; factor < factors; factor++) {
				wordFactors[word][factor] = r.nextDouble();
			}
		}

	}

	public void reComputeFactors(Corpus corpus, int factors) {
		Arrays.fill(docFactorSum, 0);
		Arrays.fill(citationFactorSum, 0);
		Arrays.fill(wordFactorSum, 0);
		for (int doc = 0; doc < corpus.docs.size(); doc++) {
			Arrays.fill(ApproximatedValue[doc], 0);
			Arrays.fill(ApproximatedWordValue[doc], 0);
			Arrays.fill(interDocFactors[doc], 0);
			Arrays.fill(interDocWordFactors[doc], 0);
		}
		for (int citation = 0; citation < corpus.docs.size(); citation++) {
			Arrays.fill(interCitationFactors[citation], 0);
		}
		for (int word = 0; word < corpus.vocabulary.size(); word++) {
			Arrays.fill(interWordFactors[word], 0);
		}

		for (int doc = 0; doc < corpus.docs.size(); doc++) {
			for (int factor = 0; factor < factors; factor++) {
				docFactorSum[factor] += documentFactors[doc][factor];
			}
		}
		for (int citation = 0; citation < corpus.docs.size(); citation++) {
			for (int factor = 0; factor < factors; factor++) {
				citationFactorSum[factor] += citationFactors[citation][factor];
			}
		}
		for (int word = 0; word < corpus.vocabulary.size(); word++) {
			for (int factor = 0; factor < factors; factor++) {
				wordFactorSum[factor] += wordFactors[word][factor];
			}
		}

		for (int doc = 0; doc < corpus.docs.size(); doc++) {
			for (int citation = 0; citation < corpus.docs.size(); citation++) {
				for (int factor = 0; factor < factors; factor++) {
					ApproximatedValue[doc][citation] += documentFactors[doc][factor]
							* citationFactors[citation][factor];
				}
			}
			for (int word = 0; word < corpus.vocabulary.size(); word++) {
				for (int factor = 0; factor < factors; factor++) {
					ApproximatedWordValue[doc][word] += documentFactors[doc][factor]
							* wordFactors[word][factor];
				}
			}
		}
		for (int doc = 0; doc < corpus.docs.size(); doc++) {
			ContextDocument document = (ContextDocument) corpus.docs.get(doc);
			int[] keys = document.citationSet.keys();
			for (int citation = 0; citation < keys.length; citation++) {
				for (int factor = 0; factor < factors; factor++) {

					double value = 1;

					/*
					 * if (corpus.citationAlphabet.contains(citation)) { value =
					 * document.citationSet .get((Integer)
					 * corpus.citationAlphabet .lookupObject(citation)); }
					 */
					if (ApproximatedValue[doc][keys[citation]] > 0) {
						interCitationFactors[keys[citation]][factor] += (documentFactors[doc][factor]
								* value / ApproximatedValue[doc][keys[citation]]);
						interDocFactors[doc][factor] += (citationFactors[keys[citation]][factor]
								* value / ApproximatedValue[doc][keys[citation]]);
					}
				}
			}
			keys = document.word_counts.keys();
			for (int word = 0; word < keys.length; word++) {
				for (int factor = 0; factor < factors; factor++) {
					if (ApproximatedWordValue[doc][keys[word]] > 0) {
						double value = document.word_counts.get(keys[word]);
						interWordFactors[keys[word]][factor] += (documentFactors[doc][factor]
								* value / ApproximatedWordValue[doc][keys[word]]);
						interDocWordFactors[doc][factor] += (wordFactors[keys[word]][factor]
								* value / ApproximatedWordValue[doc][keys[word]]);
					}
				}
			}
		}

	}

	public void learnFactors(Corpus corpus, int factors, boolean isTraining) {
		for (int doc = 0; doc < corpus.docs.size(); doc++) {
			for (int factor = 0; factor < factors; factor++) {
				if (isTraining) {
					if (citationFactorSum[factor] > 0
							&& wordFactorSum[factor] > 0)
						documentFactors[doc][factor] *= (0.9 * (interDocFactors[doc][factor] / citationFactorSum[factor]) + 0.1 * (interDocWordFactors[doc][factor] / wordFactorSum[factor]));
				} else {
					if (wordFactorSum[factor] > 0) {
						documentFactors[doc][factor] *= (interDocWordFactors[doc][factor] / wordFactorSum[factor]);
					}
				}
			}
		}
		if (isTraining) {
			for (int citation = 0; citation < corpus.docs.size(); citation++) {
				for (int factor = 0; factor < factors; factor++) {
					if (docFactorSum[factor] > 0) {
						citationFactors[citation][factor] *= (interCitationFactors[citation][factor] / docFactorSum[factor]);
					}
				}
			}
		}
		for (int word = 0; word < corpus.vocabulary.size(); word++) {
			for (int factor = 0; factor < factors; factor++) {
				if (docFactorSum[factor] > 0) {
					wordFactors[word][factor] *= (interWordFactors[word][factor] / docFactorSum[factor]);
				}
			}
		}
	}

	public double[] CalculateRecall(Corpus corpus, int K, int numtopics,
			TIntIntHashMap trainCitations) {

		double[] p_at_k = new double[K];
		// train here as well

		// initialize docs factors here
		for (int doc = 0; doc < corpus.docs.size(); doc++) {
			for (int factor = 0; factor < numtopics; factor++) {
				documentFactors[doc][factor] = r.nextDouble();
			}
		}

		for (int iter = 0; iter < 20; iter++) {
			reComputeFactors(corpus, numTopics);
			learnFactors(corpus, numTopics, false);

		}
		for (int doc = 0; doc < corpus.docs.size(); doc++) {
			ContextDocument document = (ContextDocument) corpus.docs.get(doc);
			if (document.citationSetTest.size() == 0) {
				continue;
			}
			TreeSet<DescSorter> sortedCitation = new TreeSet<DescSorter>();
			// train here on testing set as well
			int[] keys = trainCitations.keys();
			for (int citation = 0; citation < keys.length; citation++) {

				double value = 0;
				for (int t = 0; t < numtopics; t++) {// to do : change
					value += documentFactors[doc][t]
							* citationFactors[keys[citation]][t];
				}
				// System.out.println(corpus.citationAlphabet.lookupObject(
				// citation).toString()
				// + ":" + value);
				String citationid = new Integer(keys[citation]).toString();
				sortedCitation.add(new DescSorter(citationid, value));

			}
			Iterator<DescSorter> iter = sortedCitation.iterator();
			document = (ContextDocument) corpus.docs.get(doc);
			int k = 0, found = 0;
			while (iter.hasNext() && k < K) {
				DescSorter desc = iter.next();
				if (document.citationSetTest.contains(Integer.parseInt(desc
						.getID()))) {
					found++;
				}
				p_at_k[k] += ((double) found / (document.citationSetTest.size() * (corpus.docs
						.size())));
				k++;
			}
		}
		return p_at_k;
	}

	public double[] CalculatePrecision(Corpus corpus, int K, int numtopics,
			TIntIntHashMap trainCitations) {

		double[] p_at_k = new double[K];
		// train here as well

		// initialize docs factors here
		documentFactors = new double[corpus.docs.size()][numtopics];
		for (int doc = 0; doc < corpus.docs.size(); doc++) {
			for (int factor = 0; factor < numtopics; factor++) {
				documentFactors[doc][factor] = r.nextDouble();
			}
		}

		for (int iter = 0; iter < 20; iter++) {
			reComputeFactors(corpus, numTopics);
			learnFactors(corpus, numTopics, false);

		}
		for (int doc = 0; doc < corpus.docs.size(); doc++) {
			ContextDocument document = (ContextDocument) corpus.docs.get(doc);
			if (document.citationSetTest.size() == 0) {
				continue;
			}

			TreeSet<DescSorter> sortedCitation = new TreeSet<DescSorter>();
			int[] keys = trainCitations.keys();
			for (int citation = 0; citation < keys.length; citation++) {// citations
																		// =
																		// citations
																		// from
																		// test
																		// set

				double value = 0;

				for (int t = 0; t < numtopics; t++) {// to do : change
					value += documentFactors[doc][t]
							* citationFactors[keys[citation]][t];
				}

				// System.out.println(keys[citation] + ":" + value);
				String citationid = new Integer(keys[citation]).toString();
				sortedCitation.add(new DescSorter(citationid, value));

			}
			Iterator<DescSorter> iter = sortedCitation.iterator();
			document = (ContextDocument) corpus.docs.get(doc);
			int k = 0, found = 0;
			while (iter.hasNext() && k < K) {
				DescSorter desc = iter.next();

				if (document.citationSetTest.contains(Integer.parseInt(desc
						.getID()))) {
					found++;
				}
				p_at_k[k] += ((double) found / ((k + 1) * (corpus.docs.size())));
				k++;

			}
		}
		return p_at_k;
	}

	public double CalculateNDCG(Corpus corpus, int K, int numtopics,
			TIntIntHashMap trainCitations) {

		double CNDCG = 0;
		// train here as well

		// initialize docs factors here
		documentFactors = new double[corpus.docs.size()][numtopics];
		for (int doc = 0; doc < corpus.docs.size(); doc++) {
			for (int factor = 0; factor < numtopics; factor++) {
				documentFactors[doc][factor] = r.nextDouble();
			}
		}

		for (int iter = 0; iter < 20; iter++) {
			reComputeFactors(corpus, numTopics);
			learnFactors(corpus, numTopics, false);

		}
		for (int doc = 0; doc < corpus.docs.size(); doc++) {
			ContextDocument document = (ContextDocument) corpus.docs.get(doc);
			if (document.citationSetTest.size() == 0) {
				continue;
			}

			TreeSet<DescSorter> sortedCitation = new TreeSet<DescSorter>();
			int[] keys = trainCitations.keys();
			for (int citation = 0; citation < keys.length; citation++) {// citations
																		// =
																		// citations
																		// from
																		// test
																		// set

				double value = 0;

				for (int t = 0; t < numtopics; t++) {// to do : change
					value += documentFactors[doc][t]
							* citationFactors[keys[citation]][t];
				}

				// System.out.println(keys[citation] + ":" + value);
				String citationid = new Integer(keys[citation]).toString();
				sortedCitation.add(new DescSorter(citationid, value));

			}
			Iterator<DescSorter> iter = sortedCitation.iterator();
			document = (ContextDocument) corpus.docs.get(doc);
			int k = 1, found = 0;
			double prefectDCG = 1, DCG = 0;
			for (int i = 2; i <= document.citationSetTest.size(); i++) {
				prefectDCG += (Math.log(2) / Math.log(i));
			}

			while (iter.hasNext() && found < document.citationSetTest.size()) {
				DescSorter desc = iter.next();

				if (document.citationSetTest.contains(Integer.parseInt(desc
						.getID()))) {
					if (k == 1) {
						DCG += 1;
					} else {
						DCG += (Math.log(2) / Math.log(k));
					}
					found++;
				}
				k++;

			}
			CNDCG += ((double) DCG / (prefectDCG * (corpus.docs.size())));
		}
		return CNDCG;
	}

	public void sampleCorpus(Corpus corpus, int numtopics, int iterations,
			int modelType, int numFolds, int K, String output) {

		corpus.prepareSplits(numFolds, false);

		double avg_perp = 0;
		double[] p_at_k = new double[K];
		double avg_rec = 0;
		double[] r_at_k = new double[K];
		double[] NDCG = new double[numFolds];
		double avgNDCG = 0;
		for (int fold = 0; fold < 1; fold++) {
			initializeFactors(corpus, numtopics);
			corpus.docs = corpus.train_docs.get(fold);
			for (int iter = 0; iter < iterations; iter++) {
				System.out.println("Iteration=" + iter);

				reComputeFactors(corpus, numtopics);
				learnFactors(corpus, numtopics, true);
			}
			corpus.docs = corpus.test_docs.get(fold);
			corpus.isPerplex = true;
			double[] p_at_k_tmp = CalculatePrecision(corpus, K, numtopics,
					corpus.trainCitations[fold]);
			double[] r_at_k_tmp = CalculateRecall(corpus, K, numtopics,
					corpus.trainCitations[fold]);
			for (int k = 0; k < K; k++) {
				System.out.println("P@" + k + ":" + p_at_k_tmp[k]);
				System.out.println("R@" + k + ":" + r_at_k_tmp[k]);
			}
			NDCG[fold] = CalculateNDCG(corpus, K, numtopics,
					corpus.trainCitations[fold]);
			System.out.println("NDCG:" + NDCG[fold]);
			avgNDCG += NDCG[fold] / numFolds;
			for (int k = 0; k < K; k++) {
				p_at_k[k] += (p_at_k_tmp[k] / numFolds);
			}

			for (int k = 0; k < K; k++) {
				r_at_k[k] += (r_at_k_tmp[k] / numFolds);
			}

		}
		try {
			// Create file
			FileWriter fstream = new FileWriter(output);
			BufferedWriter out = new BufferedWriter(fstream);

			for (int k = 0; k < K; k++) {
				out.write("P@" + k + ":" + p_at_k[k] + "\n");
				out.write("R@" + k + ":" + r_at_k[k] + "\n");
			}
			out.write("NDCG:" + avgNDCG + "\n");
			// Close the output stream
			out.close();
		} catch (Exception e) {// Catch exception if any
			System.err.println("Error: " + e.getMessage());
		}

	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		// Read documents
		String input = args[0];
		String output = args[1];

		int windowLength = Integer.parseInt(args[2]);
		int numtopics = Integer.parseInt(args[3]);
		int iterations = Integer.parseInt(args[4]);
		int modelType = Integer.parseInt(args[5]);// 0-link, 1-cite, 2-linkplsa,
													// 3-citeplsa
		boolean isContextAware = false;
		if (modelType % 2 == 1) {
			isContextAware = true;
		}
		int numFolds = Integer.parseInt(args[6]);
		int K = Integer.parseInt(args[7]);
		double[][] plot_citing = new double[windowLength][2];
		double[][] plot_cited = new double[windowLength][2];
		double[][] plot_citing_cited = new double[windowLength][2];
		double[] citingavg = new double[windowLength];
		double[] citedavg = new double[windowLength];
		double[] citing_cited_avg = new double[windowLength];
		for (int i = 0; i < windowLength; i++) {
			citing_cited_avg[i] = 0;
		}
		NMF citationPrediction = new NMF();
		System.out.print("Assigning Ids.....");
		Corpus corpus = Corpus.assignDocumentIds(input);
		System.out.println("Done");
		System.out.print("Reading Data.....");
		corpus.readData(input);
		System.out.println("Done");
		System.out.print("Setting Window Length.....");
		corpus.setWindowLength(windowLength);
		System.out.println("Done");
		System.out.println("Size of docs=" + corpus.docs.size());
		// prepare citation alphabet

		for (int doc_index = 0; doc_index < corpus.docs.size(); doc_index++) {
			ContextDocument doc = (ContextDocument) corpus.docs.get(doc_index);
			for (int i = 0; i < doc.sentenseLength; i++) {
				if (doc.contextObject.contains(i)) {
					for (int citations = 0; citations < doc.contextObject
							.get(i).size(); citations++) {
						corpus.citationAlphabet.lookupIndex(doc.contextObject
								.get(i).get(citations), true);// initialize the
																// citation
					}
				}
			}
		}

		citationPrediction.sampleCorpus(corpus, numtopics, iterations,
				modelType, numFolds, K, output);

		// }

	}

}

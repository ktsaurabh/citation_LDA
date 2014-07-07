package edu.psu.Experiments;

import java.util.Iterator;
import java.util.Random;
import java.util.TreeSet;
import java.util.Vector;
import java.io.*;

import gnu.trove.*;

import edu.psu.topic.Model;
import edu.psu.types.*;
import edu.psu.util.*;

public class relevantCitation {

	/**
	 * @param args
	 */
	public int numTopics;
	public Random r = new Random();
	public Randoms random;
	public int rad;

	public Vector<TIntObjectHashMap<TIntDoubleHashMap>> getContext(
			ContextDocument doc) {
		Vector<TIntObjectHashMap<TIntDoubleHashMap>> all_cxt = new Vector<TIntObjectHashMap<TIntDoubleHashMap>>();
		TIntObjectHashMap<TIntDoubleHashMap> cxt_w_id = new TIntObjectHashMap<TIntDoubleHashMap>();
		int radius = rad;
		int[] words = doc.wordSequence.getFeatures();
		for (int i = 0; i < doc.sentenseLength; i++) {
			if (doc.contextObject.contains(i)) {
				TIntDoubleHashMap cxt = new TIntDoubleHashMap();
				int begin = i - radius, end = i + radius;
				// generate context vector
				if (i < radius) {
					begin = i;
				}
				if (i + radius >= doc.sentenseLength) {
					end = doc.sentenseLength - 1;
				}
				for (int j = begin; j <= end; j++) {
					for (int k = doc.SentenseBoundaries.get(j).firstElement(); k < doc.SentenseBoundaries
							.get(j).lastElement(); k++) {
						cxt.adjustOrPutValue(words[k], 1, 1);

					}
				}

				double norm = 0;

				int[] keys = cxt.keys();
				for (int j = 0; j < keys.length; j++) {
					norm += Math.pow(cxt.get(keys[j]), 2);
				}
				norm = Math.pow(norm, 0.5);

				for (int j = 0; j < keys.length; j++) {
					cxt.put(keys[j], cxt.get(keys[j]) / norm);
				}
				for (int j = 0; j < doc.contextObject.get(i).size(); j++) {
					cxt_w_id.put(doc.contextObject.get(i).get(j), cxt);// note:
																		// shallow
																		// copies
				}
				all_cxt.add(cxt_w_id);
			}
		}
		return all_cxt;
	}

	public double[] CalculateRecall(Corpus corpus, int K, int numtopics,
			TIntIntHashMap trainCitations) {
		double[] p_at_k = new double[K];
		int[] keys = trainCitations.keys();
		/*
		 * for (int i = 0; i < keys.length; i++) { ConceptVectors cv =
		 * Corpus.concepts.get(keys[i]); if (cv.contexts != null) {
		 * System.out.println(cv.docName + ":" + cv.contexts.size()); // print
		 * int[] context_keys = cv.content.keys(); for (int k = 0; k <
		 * context_keys.length; k++) { System.out.print(corpus.vocabulary
		 * .lookupObject(context_keys[k]) + " "); } System.out.println(); for
		 * (int j = 0; j < cv.contexts.size(); j++) { TIntDoubleHashMap context
		 * = cv.contexts.get(j); int[] context_keys = context.keys(); for (int k
		 * = 0; k < context_keys.length; k++) {
		 * System.out.print(corpus.vocabulary .lookupObject(context_keys[k]) +
		 * " "); } System.out.println(); } } }
		 */
		int total_seen = 0;

		for (int doc = 0; doc < corpus.docs.size(); doc++) {
			ContextDocument document = (ContextDocument) corpus.docs.get(doc);
			if (document.citationSetTest.size() == 0) {
				continue;
			}
			// Check the input here

			TreeSet<DescSorter> sortedCitation = new TreeSet<DescSorter>();
			Vector<TIntObjectHashMap<TIntDoubleHashMap>> all_cxt = getContext(document);
			for (int contexts = 0; contexts < all_cxt.size(); contexts++) {
				total_seen++;
				TIntObjectHashMap<TIntDoubleHashMap> cxt_w_id = all_cxt
						.get(contexts);
				int cited = cxt_w_id.keys()[0];
				TIntDoubleHashMap cxt = cxt_w_id.get(cited);
				keys = trainCitations.keys();
				for (int citation = 0; citation < keys.length; citation++) {
					double value = 0;
					ConceptVectors cv = Corpus.concepts.get(keys[citation]);
					// int[] content_keys =
					// Corpus.concepts.get(doc).content.keys();
					int[] content_keys = cxt.keys();
					int total_size=1;
					if (Corpus.concepts.get(keys[citation]).contexts != null && Corpus.concepts.get(keys[citation]).content!=null) {
						total_size+=Corpus.concepts.get(keys[citation]).contexts
								.size();
						for (int i = 0; i < content_keys.length; i++) {
							value += Corpus.concepts.get(keys[citation]).content
									.get(content_keys[i])
									* cxt.get(content_keys[i]);
						}
						for (int context = 0; context < Corpus.concepts
								.get(keys[citation]).contexts.size(); context++) {
							for (int i = 0; i < content_keys.length; i++) {
								value += Corpus.concepts.get(keys[citation]).contexts
										.get(context).get(content_keys[i])
										* cxt.get(content_keys[i]);
							}
						}
					}
					value /= (total_size);
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
					// if
					// (document.citationSetTest.contains(Integer.parseInt(desc.getID())))
					// {
					if (cited == Integer.parseInt(desc.getID())) {
						found++;
					}
					p_at_k[k] += ((double) found / (document.citationSetTest
							.size() * (corpus.docs.size())));
					k++;
				}
			}
		}
		for (int i = 0; i < p_at_k.length; i++) {
			p_at_k[i] /= total_seen;
		}
		return p_at_k;
	}

	public double[] CalculatePrecision(Corpus corpus, int K, int numtopics,
			TIntIntHashMap trainCitations) {
		double[] p_at_k = new double[K];
		int[] keys = trainCitations.keys();
		/*
		 * for (int i = 0; i < keys.length; i++) { ConceptVectors cv =
		 * Corpus.concepts.get(keys[i]); if (cv.contexts != null) {
		 * System.out.println(cv.docName + ":" + cv.contexts.size()); // print
		 * int[] context_keys = cv.content.keys(); for (int k = 0; k <
		 * context_keys.length; k++) { System.out.print(corpus.vocabulary
		 * .lookupObject(context_keys[k]) + " "); } System.out.println(); for
		 * (int j = 0; j < cv.contexts.size(); j++) { TIntDoubleHashMap context
		 * = cv.contexts.get(j); int[] context_keys = context.keys(); for (int k
		 * = 0; k < context_keys.length; k++) {
		 * System.out.print(corpus.vocabulary .lookupObject(context_keys[k]) +
		 * " "); } System.out.println(); } } }
		 */
		int total_seen = 0;

		for (int doc = 0; doc < corpus.docs.size(); doc++) {
			ContextDocument document = (ContextDocument) corpus.docs.get(doc);
			if (document.citationSetTest.size() == 0) {
				continue;
			}
			// Check the input here

			TreeSet<DescSorter> sortedCitation = new TreeSet<DescSorter>();
			Vector<TIntObjectHashMap<TIntDoubleHashMap>> all_cxt = getContext(document);
			for (int contexts = 0; contexts < all_cxt.size(); contexts++) {
				total_seen++;
				TIntObjectHashMap<TIntDoubleHashMap> cxt_w_id = all_cxt
						.get(contexts);
				int cited = cxt_w_id.keys()[0];
				TIntDoubleHashMap cxt = cxt_w_id.get(cited);
				keys = trainCitations.keys();
				for (int citation = 0; citation < keys.length; citation++) {
					double value = 0;
					ConceptVectors cv = Corpus.concepts.get(keys[citation]);
					// int[] content_keys =
					// Corpus.concepts.get(doc).content.keys();
					int[] content_keys = cxt.keys();
					int total_size=1;
					if (Corpus.concepts.get(keys[citation]).contexts != null && Corpus.concepts.get(keys[citation]).content!=null) {
						total_size+=Corpus.concepts.get(keys[citation]).contexts
								.size();
						for (int i = 0; i < content_keys.length; i++) {
							value += Corpus.concepts.get(keys[citation]).content
									.get(content_keys[i])
									* cxt.get(content_keys[i]);
						}
						for (int context = 0; context < Corpus.concepts
								.get(keys[citation]).contexts.size(); context++) {
							for (int i = 0; i < content_keys.length; i++) {
								value += Corpus.concepts.get(keys[citation]).contexts
										.get(context).get(content_keys[i])
										* cxt.get(content_keys[i]);
							}
						}
					}
					value /= (total_size);
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
					// if
					// (document.citationSetTest.contains(Integer.parseInt(desc.getID())))
					// {
					if (cited == Integer.parseInt(desc.getID())) {
						found++;
					}
					p_at_k[k] += ((double) found / ((k + 1)));
					k++;
				}
			}
		}
		for (int i = 0; i < p_at_k.length; i++) {
			p_at_k[i] /= total_seen;
		}
		return p_at_k;
	}

	public double CalculateNDCG(Corpus corpus, int K, int numtopics,
			TIntIntHashMap trainCitations) {
		double CNDCG = 0;
		int[] keys = trainCitations.keys();
		/*
		 * for (int i = 0; i < keys.length; i++) { ConceptVectors cv =
		 * Corpus.concepts.get(keys[i]); if (cv.contexts != null) {
		 * System.out.println(cv.docName + ":" + cv.contexts.size()); // print
		 * int[] context_keys = cv.content.keys(); for (int k = 0; k <
		 * context_keys.length; k++) { System.out.print(corpus.vocabulary
		 * .lookupObject(context_keys[k]) + " "); } System.out.println(); for
		 * (int j = 0; j < cv.contexts.size(); j++) { TIntDoubleHashMap context
		 * = cv.contexts.get(j); int[] context_keys = context.keys(); for (int k
		 * = 0; k < context_keys.length; k++) {
		 * System.out.print(corpus.vocabulary .lookupObject(context_keys[k]) +
		 * " "); } System.out.println(); } } }
		 */
		int total_seen = 0;

		for (int doc = 0; doc < corpus.docs.size(); doc++) {
			ContextDocument document = (ContextDocument) corpus.docs.get(doc);
			if (document.citationSetTest.size() == 0) {
				continue;
			}
			// Check the input here

			TreeSet<DescSorter> sortedCitation = new TreeSet<DescSorter>();
			Vector<TIntObjectHashMap<TIntDoubleHashMap>> all_cxt = getContext(document);
			for (int contexts = 0; contexts < all_cxt.size(); contexts++) {
				total_seen++;
				TIntObjectHashMap<TIntDoubleHashMap> cxt_w_id = all_cxt
						.get(contexts);
				int cited = cxt_w_id.keys()[0];
				TIntDoubleHashMap cxt = cxt_w_id.get(cited);
				keys = trainCitations.keys();
				for (int citation = 0; citation < keys.length; citation++) {
					double value = 0;
					ConceptVectors cv = Corpus.concepts.get(keys[citation]);
					// int[] content_keys =
					// Corpus.concepts.get(doc).content.keys();
					int[] content_keys = cxt.keys();
					int total_size=1;
					if (Corpus.concepts.get(keys[citation]).contexts != null && Corpus.concepts.get(keys[citation]).content!=null) {
						total_size+=Corpus.concepts.get(keys[citation]).contexts
								.size();
						for (int i = 0; i < content_keys.length; i++) {
							value += Corpus.concepts.get(keys[citation]).content
									.get(content_keys[i])
									* cxt.get(content_keys[i]);
						}
						for (int context = 0; context < Corpus.concepts
								.get(keys[citation]).contexts.size(); context++) {
							for (int i = 0; i < content_keys.length; i++) {
								value += Corpus.concepts.get(keys[citation]).contexts
										.get(context).get(content_keys[i])
										* cxt.get(content_keys[i]);
							}
						}
					}
					value /= (total_size);
					// System.out.println(corpus.citationAlphabet.lookupObject(
					// citation).toString()
					// + ":" + value);
					String citationid = new Integer(keys[citation]).toString();
					sortedCitation.add(new DescSorter(citationid, value));
				}
				Iterator<DescSorter> iter = sortedCitation.iterator();
				document = (ContextDocument) corpus.docs.get(doc);
				int k = 1;
				double DCG = 0;
				while (iter.hasNext()) {
					DescSorter desc = iter.next();
					// if
					// (document.citationSetTest.contains(Integer.parseInt(desc.getID())))
					// {
					if (cited == Integer.parseInt(desc.getID())) {
						if (k == 1) {
							DCG = 1;
						} else {
							DCG = (Math.log(2) / Math.log(k));
						}
						break;
					}

					k++;
				}

				CNDCG += DCG;
			}
			CNDCG /= total_seen;
		}

		return CNDCG;
	}

	public void sampleCorpus(Corpus corpus, int numtopics, int iterations,
			int modelType, int numFolds, int K, String output) {

		corpus.prepareSplits(numFolds, true);
		double[] p_at_k = new double[K];
		double[] r_at_k = new double[K];
		double[] NDCG = new double[numFolds];
		double avgNDCG = 0;
		for (int fold = 0; fold < numFolds; fold++) {

			corpus.docs = corpus.train_docs.get(fold);

			/*
			 * for (int iter = 0; iter < iterations; iter++) { int sampleSize =
			 * 0;
			 * 
			 * System.out.println("Samples=" + sampleSize);
			 * System.out.println("Iter=" + iter); //
			 * model.estimateParameters(corpus); // if (iter > 10 && iter % 3 ==
			 * 0) { // model.estimateParameters(corpus); // } //
			 * System.out.println("Loglikelihood=" // + empiricalLikelihood(10,
			 * corpus)); }
			 */
			corpus.docs = corpus.test_docs.get(fold);
			corpus.isPerplex = true;

			double perp_report = Double.MAX_VALUE;

			for (int iter = iterations; iter < iterations + 30; iter++) {

			}
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
		if (args.length < 7) {
			System.out
					.println("Not enough parameters; Please see readme file for parameter values.");
			System.exit(0);
		}
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

		double[] citing_cited_avg = new double[windowLength];
		for (int i = 0; i < windowLength; i++) {
			citing_cited_avg[i] = 0;
		}
		relevantCitation relevantcitation = new relevantCitation();
		relevantcitation.rad = windowLength;
		System.out.print("Assigning Ids.....");
		Corpus corpus = Corpus.assignDocumentIds(input);
		System.out.println("Done");
		System.out.print("Reading Data.....");
		corpus.readData(input);
		/*
		 * TreeSet<IDSorter> sortedCitation = new TreeSet<IDSorter>(); for(int
		 * i=0;i<corpus.docs.size();i++){ ContextDocument doc =
		 * (ContextDocument) corpus.docs.get(i); sortedCitation.add(new
		 * IDSorter(i,(double)doc.citationSet.size()));
		 * 
		 * } Iterator<IDSorter> iter = sortedCitation.iterator();
		 * 
		 * 
		 * while (iter.hasNext()) { IDSorter desc = iter.next();
		 * System.out.println(desc.getID() + ":" + desc.getWeight()); }
		 * System.exit(0);
		 */

		System.out.println("Done");
		int count = 0;
		/*
		 * for (int i = 0; i < Corpus.concepts.size(); i++) { ConceptVectors cv
		 * = Corpus.concepts.get(i); if (cv.contexts != null) { count++;
		 * System.out.println(cv.docName + ":" + cv.contexts.size()); //print
		 * contexts for(int j=0;j<cv.contexts.size();j++){ TIntDoubleHashMap
		 * context = cv.contexts.get(j); int[] keys= context.keys(); for(int
		 * k=0;k<keys.length;k++){
		 * System.out.print(corpus.vocabulary.lookupObject(keys[k])+" "); }
		 * System.out.println(); } } }
		 */
		// System.out.print(count);
		// System.out.print("Setting Window Length.....");
		/*
		 * corpus.setWindowLength(windowLength); System.out.println("Done");
		 * System.out.println("Size of docs=" + corpus.docs.size());
		 */
		// System.exit(0);
		relevantcitation.sampleCorpus(corpus, numtopics, iterations, modelType,
				numFolds, K, output);

		System.exit(0);

	}

}

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

public class CitationPrediction {

	/**
	 * @param args
	 */
	public int numTopics;
	public Random r = new Random();
	public Randoms random;
	public boolean ifAdaptive;

	public double[] CalculateRecall(Corpus corpus, int K, int numtopics,
			TIntIntHashMap trainCitations) {

		double[] p_at_k = new double[K];
		for (int doc = 0; doc < corpus.docs.size(); doc++) {
			ContextDocument document = (ContextDocument) corpus.docs.get(doc);
			if (document.citationSetTest.size() == 0) {
				continue;
			}
			int[] keys = trainCitations.keys();
			TreeSet<DescSorter> sortedCitation = new TreeSet<DescSorter>();
			for (int citation = 0; citation < keys.length; citation++) {

				double value = 0;
				for (int t = 0; t < numtopics; t++) {// to do : change
					value += corpus.theta_train[doc][t]
							* corpus.psi_train[keys[citation]][t];
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
			while (iter.hasNext() && k < K && found < document.citationSetTest.size()) {
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
		System.out.println("Stat: test docs:" + corpus.docs.size()
				+ ":Total citations in training:" + trainCitations.size());
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
					value += corpus.theta_train[doc][t]
							* corpus.psi_train[keys[citation]][t];
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
			while (iter.hasNext() && k < K && found < document.citationSetTest.size()) {
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
		System.out.println("Stat: test docs:" + corpus.docs.size()
				+ ":Total citations in training:" + trainCitations.size());
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
					value += corpus.theta_train[doc][t]
							* corpus.psi_train[keys[citation]][t];
				}
				// System.out.println(corpus.citationAlphabet.lookupObject(
				// citation).toString()
				// + ":" + value);
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
		Model model = null;
		boolean isContextAware = false;
		if (modelType % 2 == 1) {
			isContextAware = true;
		}
		switch (modelType) {
		case 0:
			model = Model.factory("linklda");
			break;
		case 1:
			model = Model.factory("citelda");
			break;
		case 2:
			model = Model.factory("linkplsalda");
			break;
		case 3:
			model = Model.factory("citeplsalda");
			break;
		default:
			System.out.println("Enter number between 0-3");
			System.exit(0);

		}

		model.numTopics = numtopics;
		corpus.prepareSplits(numFolds, false);
		Vector<Document> docs = corpus.docs;
		double[] p_at_k = new double[K];
		double[] r_at_k = new double[K];
		double[] NDCG = new double[numFolds];
		double avgNDCG = 0;
		for (int fold = 0; fold < numFolds; fold++) {
			model.numSamples = 0;
			corpus.docs = docs;
			model.InitializeParameters(corpus);
			corpus.docs = corpus.train_docs.get(fold);
			model.InitializeAssignments(corpus, new LabelAlphabet());
			for (int iter = 0; iter < iterations; iter++) {
				int sampleSize = 0;
				model.sampleCorpus(corpus, iter, isContextAware);
				if (iter > 5 && iter < 15 && iter % 2 == 0 && ifAdaptive
						&& fold == 0)
					corpus.adaptiveWindowLength(fold);
				// System.out.println("Samples=" + sampleSize);
				System.out.println("Iter=" + iter);
				if(iter>0 && iter%20==0){
					model.learnParameters(corpus);
					//for(int i=0;i<model.numTopics;i++){
						//System.out.print(corpus.alpha[i]+"\t");
					//}
				}
				// model.estimateParameters(corpus);
				// if (iter > 10 && iter % 3 == 0) {
				// model.estimateParameters(corpus);
				// }
				// System.out.println("Loglikelihood="
				// + empiricalLikelihood(10, corpus));
			}
			//printTopicDiagnostics(corpus);
			corpus.docs = corpus.test_docs.get(fold);
			corpus.isPerplex = true;
			model.InitializeAssignments(corpus, new LabelAlphabet());
			double perp_report = Double.MAX_VALUE;

			for (int iter = iterations; iter < iterations + 30; iter++) {
				model.estimateParameters(corpus);
				for (int i = 0; i < corpus.docs.size(); i++) {
					model.sampleOneDocument(corpus,
							(ContextDocument) corpus.docs.get(i));
				}

				// link prediction

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
	private class Kindall_Statistics {
		double max, min, mean, variance, min_index, max_index;
	}

	public Kindall_Statistics compute_kindall_statistics(Vector<Double> kindalls) {
		Kindall_Statistics kindall_statistics = new Kindall_Statistics();

		double mean = 0, max = Double.MIN_VALUE, min = Double.MAX_VALUE;
		for (int i = 0; i < kindalls.size(); i++) {
			if (max < kindalls.get(i)) {
				max = kindalls.get(i);
				kindall_statistics.max_index = i;
			}
			if (min > kindalls.get(i)) {
				min = kindalls.get(i);
				kindall_statistics.min_index = i;
			}
			mean += kindalls.get(i);
		}
		mean /= kindalls.size();
		kindall_statistics.mean = mean;
		kindall_statistics.min = min;
		kindall_statistics.max = max;

		double variance = 0;
		for (int i = 0; i < kindalls.size(); i++) {
			variance += Math.pow(mean - kindalls.get(i), 2);
		}
		variance = Math.pow(variance, 0.5);
		variance /= kindalls.size();
		kindall_statistics.variance = variance;

		return kindall_statistics;
	}

	public void printTopicDiagnostics(Corpus corpus) {
		System.out.println("Printing Topic Kindall-Taus.....");
		Vector<RankedList> cited = new Vector<RankedList>();
		
		Vector<RankedList> topics = new Vector<RankedList>();
		RankedList kindall = new RankedList();
		for (int i = 0; i < numTopics; i++) {
			Vector<Integer> top_100_cited = new Vector<Integer>();
			Vector<Integer> top_100_citing = new Vector<Integer>();
			RankedList words = new RankedList();
			RankedList citedAuthors = new RankedList();
			
			for (int k = 0; k < corpus.vocabulary.size(); k++) {
				words.add(k, corpus.phi_train[k][i]);
			}
			words.sortListDesc();
			topics.add(words);
			// collect top 100 authors in citing and cited set
			int collected = 0;
			TreeSet<IDSorter> sortedauthors = new TreeSet<IDSorter>();
			for (int k = 0; k < corpus.docAlphabet.size(); k++) {
				if (!corpus.docAlphabet.lookupObject(k).equals("dummy")) {
					sortedauthors.add(new IDSorter(k, corpus.psi_train[k][i]));
				}
			}
			Iterator<IDSorter> iter = sortedauthors.iterator();

			while (iter.hasNext() && collected < 100) {
				IDSorter auth = iter.next();
				top_100_cited.add(auth.getID());
				if (!citedAuthors.values.containsKey(auth.getID())) {
					citedAuthors.add(auth.getID(),
							corpus.psi_train[auth.getID()][i]);
				}
				
				collected++;
			}

			
			
			// collection end

			citedAuthors.sortListDesc();
			cited.add(citedAuthors);
			

		}
		kindall.sortListDesc();

		try {
			FileWriter fw = new FileWriter("./doc_topic_kindall_tau_"
					+ numTopics);
			BufferedWriter bw = new BufferedWriter(fw);
			// get global statistics for kindall taus
			Vector<Double> kindalls = new Vector<Double>();
			for (int i = 0; i < numTopics; i++) {
				kindalls.add(0.0);
			}

			
			// influential and cited authors
			for (int i = 0; i < numTopics; i++) {
				RankedList citedAuthors = cited.get(i);
				kindalls.setElementAt(citedAuthors.KindallTauCoefficient(Corpus.citedList),i);
			}
			Kindall_Statistics ks = compute_kindall_statistics(kindalls);
			bw.write("-------------------------influential and cited authors------------------------------\n");
			bw.write("Mean = " + ks.mean + ": Variance = " + ks.variance + "\n");
			bw.write("Max = " + ks.max + " for topic = " + ks.max_index + "\n");
			bw.write("Min = " + ks.min + " for topic = " + ks.min_index
					+ "\n\n");

			// interesting and citing authors

			

			bw.write("-------------------------------------Topics------------------------------------------------\n");
			for (int i = 0; i < numTopics; i++) {
				int topic = kindall.indices.get(i);
				RankedList words = topics.get(topic);
				RankedList citedAuthors = cited.get(topic);
				
				bw.write("Topic-"
						+ topic
						+ ":Kindall-tau="
						+ kindall.values.get(topic)
						+ ": Kindall Tau Cited:"
						+ citedAuthors
								.KindallTauCoefficient(Corpus.citedList)
						
						
						+ "\n");
				bw.write("---------------------------------------------------------------------------------------------------------\n");
				bw.write("Top Words\tTop Influential Papers\n");
				bw.write("---------------------------------------------------------------------------------------------------------\n");
				for (int j = 0; j < 20; j++) {
					double value = words.values.get(words.indices.get(j));
					bw.write(Corpus.vocabulary.lookupObject(words.indices
							.get(j)) + ":" + value + "\t");
					
					value = citedAuthors.values
							.get(citedAuthors.indices.get(j));
					bw.write(Corpus.docAlphabet
							.lookupObject(citedAuthors.indices.get(j))
							+ ":"
							+ value + "\n");
				}
				bw.write("---------------------------------------------------------------------------------------------------------\n");

			}
			bw.close();
			fw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
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

		CitationPrediction citationPrediction = new CitationPrediction();
		int adapt = Integer.parseInt(args[8]);
		if (adapt == 1) {
			citationPrediction.ifAdaptive = true;
		} else {
			citationPrediction.ifAdaptive = false;
		}

		double[] citing_cited_avg = new double[windowLength];
		for (int i = 0; i < windowLength; i++) {
			citing_cited_avg[i] = 0;
		}

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
		System.out.print("Setting Window Length.....");
		corpus.setWindowLength(windowLength);
		System.out.println("Done");
		System.out.println("Size of docs=" + corpus.docs.size());

		citationPrediction.sampleCorpus(corpus, numtopics, iterations,
				modelType, numFolds, K, output);

		System.exit(0);

	}

}

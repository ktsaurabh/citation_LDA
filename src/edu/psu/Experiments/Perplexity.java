package edu.psu.Experiments;

import java.util.Random;
import java.io.*;

import gnu.trove.*;

import edu.psu.topic.Model;
import edu.psu.types.*;
import edu.psu.util.*;


public class Perplexity {

	/**
	 * @param args
	 */
	public int numTopics;
	public Random r = new Random();
	public Randoms random;

	public void sampleCorpus(Corpus corpus, int numtopics, int iterations,
			int modelType, int numFolds, int K, String output, boolean ifAdaptive) {
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
		corpus.prepareSplits(numFolds, true);
		model.InitializeParameters(corpus);
		double avg_perp = 0;

		for (int fold = 0; fold < numFolds; fold++) {
			model.numSamples = 0;
			corpus.docs = corpus.train_docs.get(fold);
			model.InitializeAssignments(corpus, new LabelAlphabet());
			for (int iter = 0; iter < iterations; iter++) {
				int sampleSize = 0;
				model.sampleCorpus(corpus, iter, isContextAware);
				if (iter > 5 && iter < 15 && iter % 2 == 0 && ifAdaptive
						&& fold == 0)
					corpus.adaptiveWindowLength(fold);
				
				System.out.print("<" + iter + ">, ");
				// model.estimateParameters(corpus);
				if(iter>5 && iter%2==0)
					model.learnParameters(corpus);
				// System.out.println("Loglikelihood="
				// + empiricalLikelihood(10, corpus));
			}
			corpus.docs = corpus.test_docs.get(fold);
			corpus.isPerplex = true;
			model.InitializeAssignments(corpus, new LabelAlphabet());
			double perp_report = Double.MAX_VALUE;

			for (int iter = iterations; iter < iterations + 20; iter++) {
				model.estimateParameters(corpus);
				for (int i = 0; i < corpus.docs.size(); i++) {
					model.sampleOneDocument(corpus,
							(ContextDocument) corpus.docs.get(i));
				}
				if (iter > iterations + 10) {
					model.estimateParameters(corpus);

					double perp = model.modelLogLikelihood(corpus);// sampleLikelihood(10,
																	// corpus);//model.testPerplexity(10,
																	// corpus);
					//double perp = model.testPerplexity(10, corpus);
					//if (perp_report > perp) {
					//	perp_report = perp;
					//}
					avg_perp+=perp;
					switch (modelType) {
					case 0:
						System.out.println("Link-LDA:Perp=" + perp);
						break;
					case 1:
						System.out.println("Cink-LDA:Perp=" + perp);
						break;
					case 2:
						System.out.println("Link-PLSA-LDA:Perp=" + perp);
						break;
					case 3:
						System.out.println("Cite-PLSA-LDA:Perp=" + perp);
						break;
					default:
						System.out.println("Enter number between 0-3");
						System.exit(0);
					
					}
				}
				// link prediction

			}

		}
		avg_perp/=numFolds;
		try {
			// Create file
			FileWriter fstream = new FileWriter(output);
			BufferedWriter out = new BufferedWriter(fstream);
			out.write("Average Perplexity:" + avg_perp + "\n");

			// Close the output stream
			out.close();
		} catch (Exception e) {// Catch exception if any
			System.err.println("Error: " + e.getMessage());
		}

	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		// Read documents
		if(args.length<6){
			System.out.println("Not enough parameters; Please see readme file for parameter values.");
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
		int K = 1;//Integer.parseInt(args[7]);
		int isAdapt = Integer.parseInt(args[7]);
		boolean ifAdaptive = true;
		if(isAdapt==0){
			ifAdaptive=false;
		}
			
		
		double[] citing_cited_avg = new double[windowLength];
		for (int i = 0; i < windowLength; i++) {
			citing_cited_avg[i] = 0;
		}
		Perplexity perplexity = new Perplexity();
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

		perplexity.sampleCorpus(corpus, numtopics, iterations, modelType,
				numFolds, K, output, ifAdaptive);

		System.exit(0);
		
	}

}

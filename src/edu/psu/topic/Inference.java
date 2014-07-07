package edu.psu.topic;

import edu.psu.types.Corpus;
import edu.psu.types.Document;
import gnu.trove.TIntIntHashMap;
import gnu.trove.TIntObjectHashMap;

import java.util.Random;
import java.util.Vector;

public abstract class Inference {
	static long seed;
	static int numTopics;
	static Random r = new Random(seed);
	Vector<Document> docs;
	TIntIntHashMap global_topic_counts;// topic counts
	TIntObjectHashMap<TIntIntHashMap> word_topic_counts;// word-topic counts
	TIntObjectHashMap<TIntIntHashMap> link_topic_counts;// link-topic counts
	abstract void initializeCorpus(Corpus corpus);

}

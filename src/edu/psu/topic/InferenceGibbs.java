package edu.psu.topic;

import edu.psu.types.*;
import gnu.trove.TIntIntHashMap;
import gnu.trove.TIntObjectHashMap;

import java.util.*;

public class InferenceGibbs extends Inference {

	void initializeCorpus(Corpus corpus) {
		Vector<Document> docs = corpus.docs;
		for (Document doc : docs) {
			global_topic_counts = new TIntIntHashMap();
			word_topic_counts = new TIntObjectHashMap<TIntIntHashMap>();
			link_topic_counts = new TIntObjectHashMap<TIntIntHashMap>();
			doc.topicAssignments = new LabelSequence(Corpus.vocabulary);
			// wordSequence = new FeatureSequence(Corpus.vocabulary);
			doc.topic_counts = new TIntIntHashMap();
			// links = new Vector<Integer>();
			doc.linksTopics = new Vector<Integer>();
			// String[] words = content.split(" ");
			// String[] edges = citations.split(" ");

			int[] words = doc.wordSequence.getFeatures();
			// int[] topics = doc.topicAssignments.getFeatures();
			

			for (int i = 0; i < doc.wordSequence.size(); i++) {

				int topic = r.nextInt(numTopics);
				doc.topicAssignments.add(topic);
				doc.topic_counts.adjustOrPutValue(topic, 1, 1);
				if (word_topic_counts.containsKey(words[i])) {
					word_topic_counts.get(words[i]).adjustOrPutValue(topic, 1,
							1);
				} else {
					TIntIntHashMap topicMap = new TIntIntHashMap();
					topicMap.adjustOrPutValue(topic, 1, 1);
					word_topic_counts.put(words[i], topicMap);
				}
				global_topic_counts.adjustOrPutValue(topic, 1, 1);

			}
			for (int i = 0; i < doc.links.size(); i++) {

				int topic = r.nextInt(numTopics);
				doc.linksTopics.add(topic);
				doc.topic_counts.adjustOrPutValue(topic, 1, 1);
				if (link_topic_counts.containsKey(doc.links.get(i))) {
					link_topic_counts.get(doc.links.get(i)).adjustOrPutValue(
							topic, 1, 1);
				} else {
					TIntIntHashMap topicMap = new TIntIntHashMap();
					topicMap.adjustOrPutValue(topic, 1, 1);
					link_topic_counts.put(doc.links.get(i), topicMap);
				}
				global_topic_counts.adjustOrPutValue(topic, 1, 1);
			}
		}

	}
}

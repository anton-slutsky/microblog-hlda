package com.zelantsoft.micro_hlda;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import cc.mallet.pipe.CharSequence2TokenSequence;
import cc.mallet.pipe.CharSequenceLowercase;
import cc.mallet.pipe.Pipe;
import cc.mallet.pipe.SerialPipes;
import cc.mallet.pipe.TokenSequence2FeatureSequence;
import cc.mallet.types.FeatureSequence;
import cc.mallet.types.Instance;
import cc.mallet.types.InstanceList;
import cc.mallet.util.Randoms;


public class TopicModel {

	
	
	public static Map<Integer, Double> levelSpecAll = new HashMap<Integer, Double>();
	public static Map<Integer, Double> levelExpRankAll = new HashMap<Integer, Double>();
	public static Map<Integer, Integer> levelSpecCnt = new HashMap<Integer, Integer>();
	
	public static double[] dfs; 
	
	
	public static void main(String[] argz) throws Exception {
	
		double [] normVect;
		double magNorm;
		
		
		String dataFile = "/Users/anton/Dev/Data/tw_trec_clean_1000.txt";
		

		// Begin by importing documents from text to feature sequences
		ArrayList<Pipe> pipeList = new ArrayList<Pipe>();

		// Pipes: lowercase, tokenize, remove stopwords, map to features
		pipeList.add(new CharSequenceLowercase());
		pipeList.add(new CharSequence2TokenSequence(Pattern.compile("\\S+")));

		pipeList.add(new TokenSequence2FeatureSequence());

		final InstanceList instances = new InstanceList(new SerialPipes(
				pipeList));
		
		Reader fileReader = new InputStreamReader(new FileInputStream(new File(
				dataFile)), "UTF-8");
		instances.addThruPipe(new CsvIterator(fileReader, Pattern
				.compile("\\S+"), 1, 0, 0)); // data, label, name fields


		dfs = new double[instances.getAlphabet().size()];

		double sum = 0;
		for (Instance instance : instances) {
			FeatureSequence fs = (FeatureSequence) instance.getData();

			int seqLen = fs.getLength();
			int type, token;
			Set<Integer> seen = new HashSet<Integer>();
			for (token = 0; token < seqLen; token++) {
				type = fs.getIndexAtPosition(token);
				if (!seen.contains(type)) {
					dfs[type]++;
					sum ++;
					seen.add(type);
				}
			}
		}
		
		normVect = new double[dfs.length];
		
		double chck = 0;
		for(int i = 0; i < dfs.length; ++i)
		{
			normVect[i] = (dfs[i] + 0.05)/(sum + (0.05*dfs.length));
			chck += normVect[i]; 
		}
		
		magNorm = Math.sqrt(dot(normVect, normVect));
				
		System.out.println(chck);
		
		class C implements Comparable<C>
		{

			final int type;
			final double freq;
			
			public C(int type, double freq) {
				super();
				this.type = type;
				this.freq = freq;
			}

			@Override
			public int compareTo(C c) {
				return new Double(c.freq).compareTo(freq);
			}
			public String toString()
			{
				return "["+type+"] "+freq;
			}
		}
		List<C> list = new ArrayList<C>(dfs.length);
		for(int i = 0; i < dfs.length; ++i)
		{
			list.add(new C(i, dfs[i]));
		}
		Collections.sort(list);
		System.out.println(list);
		for(int i = 0; i < list.size(); ++i)
		{
			dfs[list.get(i).type] = i;
		}

		 MicroblogHLDA hlda = new MicroblogHLDA(dfs, normVect);
		 hlda.initialize(instances, null, 5, new Randoms());
		 hlda.estimate(250);

		
	}
	public static double dot(double [] v1, double [] v2)
    {
    	double res = 0;
    	for(int i = 0; i < v1.length; ++i)
    	{
    		res += (v1[i] * v2[i]);
    	}
    	return res;
    }

}

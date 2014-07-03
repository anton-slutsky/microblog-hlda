package com.zelantsoft.micro_hlda;

import gnu.trove.TIntIntHashMap;
import gnu.trove.TObjectDoubleHashMap;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;

import cc.mallet.types.Alphabet;
import cc.mallet.types.Dirichlet;
import cc.mallet.types.FeatureSequence;
import cc.mallet.types.IDSorter;
import cc.mallet.types.Instance;
import cc.mallet.types.InstanceList;
import cc.mallet.util.Randoms;


public class MicroblogHLDA {

	final double [] normVect;
	
    InstanceList instances;
    InstanceList testing;

    NCRPNode rootNode, node;

    int numLevels;
    int numDocuments;
    int numTypes;

    double alpha; // smoothing on topic distributions
    double gamma; // "imaginary" customers at the next, as yet unused table
    double eta;   // smoothing on word distributions
    double etaSum[];
    double base;
    double phi, phiSum;

    int[][] levels; // indexed < doc, token >
    NCRPNode[] documentLeaves; // currently selected path (ie leaf node) through the NCRP tree

	int totalNodes = 0;

	String stateFile = "hlda.state";

    Randoms random;

	boolean showProgress = true;
	
	int displayTopicsInterval = 50;
	int numWordsToDisplay = 10;

	double [] dfs;
	
    public MicroblogHLDA (double [] dfs, double [] normVect) {
		alpha = 0.1;
		gamma = 1.0;
		phi = 1.0;
		eta = 100;//0.01;
		this.dfs = dfs; 
		this.normVect = normVect;
		
    }

	public void setAlpha(double alpha) {
		this.alpha = alpha;
	}

	public void setGamma(double gamma) {
		this.gamma = gamma;
	}

	public void setEta(double eta) {
		this.eta = eta;
	}

	public void setStateFile(String stateFile) {
		this.stateFile = stateFile;
	}

	public void setTopicDisplay(int interval, int words) {
		displayTopicsInterval = interval;
		numWordsToDisplay = words;
	}

	/**  
	 *  This parameter determines whether the sampler outputs 
	 *   shows progress by outputting a character after every iteration.
	 */
	public void setProgressDisplay(boolean showProgress) {
		this.showProgress = showProgress;
	}

    public void initialize(InstanceList instances, InstanceList testing,
						   int numLevels, Randoms random) {
    	
    	
		this.instances = instances;
		this.testing = testing;
		this.numLevels = numLevels;
		this.random = random;
		
		double max = 0, c = 0;
		for(int i = 0; i < dfs.length; ++i)
		{
			if(dfs[i]>max)
			{
				max = dfs[i];
				c = i;
			}
		}
		System.out.println("Max is "+max+" at "+c+" of "+dfs.length);

		base = Math.pow(dfs.length+1, (1d/(numLevels)));
		
		if (! (instances.get(0).getData() instanceof FeatureSequence)) {
			throw new IllegalArgumentException("Input must be a FeatureSequence, using the --feature-sequence option when impoting data, for example");
		}

		numDocuments = instances.size();
		numTypes = instances.getDataAlphabet().size();
	
		phiSum = phi * numTypes;
	
		etaSum = new double[numLevels];
		for(int level = 0; level < numLevels; ++level)
		{
			for(int type = 0; type < instances.getAlphabet().size(); type++)
			{
				etaSum[level] += getEta(type, level);
			}
		}

		// Initialize a single path

		NCRPNode[] path = new NCRPNode[numLevels];

		rootNode = new NCRPNode(numTypes);

		levels = new int[numDocuments][];
		documentLeaves = new NCRPNode[numDocuments];

		// Initialize and fill the topic pointer arrays for 
		//  every document. Set everything to the single path that 
		//  we added earlier.
		for (int doc=0; doc < numDocuments; doc++) {
            FeatureSequence fs = (FeatureSequence) instances.get(doc).getData();
            int seqLen = fs.getLength();

			path[0] = rootNode;
			rootNode.customers++;
			for (int level = 1; level < numLevels; level++) {
				path[level] = path[level-1].select();
				path[level].customers++;
			}
			node = path[numLevels - 1];
	    
			levels[doc] = new int[seqLen];
			documentLeaves[doc] = node;

			for (int token=0; token < seqLen; token++) {
				int type = fs.getIndexAtPosition(token);
				levels[doc][token] = random.nextInt(numLevels);
				node = path[ levels[doc][token] ];
				node.totalTokens++;
				node.typeCounts[type]++;
			}
		}
		
		
	}

	public void estimate(int numIterations) {
		for (int iteration = 1; iteration <= numIterations; iteration++) {
			for (int doc=0; doc < numDocuments; doc++) {
				samplePath(doc, iteration);
			}
			for (int doc=0; doc < numDocuments; doc++) {
				sampleTopics(doc);
			}
			
			if (showProgress) {
				System.out.print(".");
				if (iteration % 50 == 0) {
					System.out.println(" " + iteration);
				}
			}

			if (iteration % displayTopicsInterval == 0) {
				printNodes();
			}
		}
    }

    public void samplePath(int doc, int iteration) {
		NCRPNode[] path = new NCRPNode[numLevels];
		NCRPNode node;
		int level, token, type;

		node = documentLeaves[doc];
		for (level = numLevels - 1; level >= 0; level--) {
			path[level] = node;
			node = node.parent;
		}

		documentLeaves[doc].dropPath();

		TObjectDoubleHashMap<NCRPNode> nodeWeights = 
			new TObjectDoubleHashMap<NCRPNode>();
	
		// Calculate p(c_m | c_{-m})
		calculateNCRP(nodeWeights, rootNode, 0.0);

		// Add weights for p(w_m | c, w_{-m}, z)
	
		// The path may have no further customers and therefore
		//  be unavailable, but it should still exist since we haven't
		//  reset documentLeaves[doc] yet...
	
		TIntIntHashMap[] typeCounts = new TIntIntHashMap[numLevels];

		int[] docLevels;

		for (level = 0; level < numLevels; level++) {
			typeCounts[level] = new TIntIntHashMap();
		}

		docLevels = levels[doc];
		FeatureSequence fs = (FeatureSequence) instances.get(doc).getData();
	    
		// Save the counts of every word at each level, and remove
		//  counts from the current path

		for (token = 0; token < docLevels.length; token++) {
			level = docLevels[token];
			type = fs.getIndexAtPosition(token);
	    
			if (! typeCounts[level].containsKey(type)) {
				typeCounts[level].put(type, 1);
			}
			else {
				typeCounts[level].increment(type);
			}

			path[level].typeCounts[type]--;
			assert(path[level].typeCounts[type] >= 0);
	    
			path[level].totalTokens--;	    
			assert(path[level].totalTokens >= 0);
		}

		// Calculate the weight for a new path at a given level.
		double[] newTopicWeights = new double[numLevels];
		for (level = 1; level < numLevels; level++) {  // Skip the root...
			int[] types = typeCounts[level].keys();
			int totalTokens = 0;

			for (int t: types) {
				for (int i=0; i<typeCounts[level].get(t); i++) {
					newTopicWeights[level] += 
						Math.log((getEta(t, level) + i) / (getEtaSum(t, level) + totalTokens));
					totalTokens++;
				}
			}

			//if (iteration > 1) { System.out.println(newTopicWeights[level]); }
		}
	
		calculateWordLikelihood(nodeWeights, rootNode, 0.0, typeCounts, newTopicWeights, 0, iteration);

		NCRPNode[] nodes = nodeWeights.keys(new NCRPNode[] {});
		double[] weights = new double[nodes.length];
		double sum = 0.0;
		double max = Double.NEGATIVE_INFINITY;

		// To avoid underflow, we're using log weights and normalizing the node weights so that 
		//  the largest weight is always 1.
		for (int i=0; i<nodes.length; i++) {
			if (nodeWeights.get(nodes[i]) > max) {
				max = nodeWeights.get(nodes[i]);
			}
		}

		for (int i=0; i<nodes.length; i++) {
			weights[i] = Math.exp(nodeWeights.get(nodes[i]) - max);

			/*
			  if (iteration > 1) {
			  if (nodes[i] == documentLeaves[doc]) {
			  System.out.print("* ");
			  }
			  System.out.println(((NCRPNode) nodes[i]).level + "\t" + weights[i] + 
			  "\t" + nodeWeights.get(nodes[i]));
			  }
			*/

			sum += weights[i];
		}

		//if (iteration > 1) {System.out.println();}

		node = nodes[ random.nextDiscrete(weights, sum) ];

		// If we have picked an internal node, we need to 
		//  add a new path.
		if (! node.isLeaf()) {
			node = node.getNewLeaf();
		}
	
		node.addPath();
		documentLeaves[doc] = node;

		for (level = numLevels - 1; level >= 0; level--) {
			int[] types = typeCounts[level].keys();

			for (int t: types) {
				node.typeCounts[t] += typeCounts[level].get(t);
				node.totalTokens += typeCounts[level].get(t);
			}

			node = node.parent;
		}
    }

    public void calculateNCRP(TObjectDoubleHashMap<NCRPNode> nodeWeights, 
							  NCRPNode node, double weight) {
		for (NCRPNode child: node.children) {
			calculateNCRP(nodeWeights, child,
						  weight + Math.log((double) child.customers / (node.customers + gamma)));
		}

		nodeWeights.put(node, weight + Math.log(gamma / (node.customers + gamma)));
    }

    public void calculateWordLikelihood(TObjectDoubleHashMap<NCRPNode> nodeWeights,
										NCRPNode node, double weight, 
										TIntIntHashMap[] typeCounts, double[] newTopicWeights,
										int level, int iteration) {
	
		// First calculate the likelihood of the words at this level, given
		//  this topic.
		double nodeWeight = 0.0;
		int[] types = typeCounts[level].keys();
		int totalTokens = 0;
	
		//if (iteration > 1) { System.out.println(level + " " + nodeWeight); }

		for (int type: types) {
			for (int i=0; i<typeCounts[level].get(type); i++) {
				nodeWeight +=
						Math.log((phi + node.typeCounts[type] + i) /
								 (phiSum + node.totalTokens + totalTokens));
				totalTokens++;

				/*
				  if (iteration > 1) {
				  System.out.println("(" +eta + " + " + node.typeCounts[type] + " + " + i + ") /" + 
				  "(" + etaSum + " + " + node.totalTokens + " + " + totalTokens + ")" + 
				  " : " + nodeWeight);
				  }
				*/

			}
		}

		// Propagate that weight to the child nodes
		for (NCRPNode child: node.children) {
            calculateWordLikelihood(nodeWeights, child, weight + nodeWeight,
									typeCounts, newTopicWeights, level + 1, iteration);
        }

		// Finally, if this is an internal node, add the weight of
		//  a new path

		level++;
		while (level < numLevels) {
			nodeWeight += newTopicWeights[level];
			level++;
		}

		nodeWeights.adjustValue(node, nodeWeight);

    }

    /** Propagate a topic weight to a node and all its children.
		weight is assumed to be a log.
	*/
    public void propagateTopicWeight(TObjectDoubleHashMap<NCRPNode> nodeWeights,
									 NCRPNode node, double weight) {
		if (! nodeWeights.containsKey(node)) {
			// calculating the NCRP prior proceeds from the
			//  root down (ie following child links),
			//  but adding the word-topic weights comes from
			//  the bottom up, following parent links and then 
			//  child links. It's possible that the leaf node may have
			//  been removed just prior to this round, so the current
			//  node may not have an NCRP weight. If so, it's not 
			//  going to be sampled anyway, so ditch it.
			return;
		}
	
		for (NCRPNode child: node.children) {
			propagateTopicWeight(nodeWeights, child, weight);
		}

		nodeWeights.adjustValue(node, weight);
    }

    public void sampleTopics(int doc) {
		FeatureSequence fs = (FeatureSequence) instances.get(doc).getData();
		int seqLen = fs.getLength();
		int[] docLevels = levels[doc];
		NCRPNode[] path = new NCRPNode[numLevels];
		NCRPNode node;
		int[] levelCounts = new int[numLevels];
		int type, token, level;
		double sum;

		// Get the leaf
		node = documentLeaves[doc];
		for (level = numLevels - 1; level >= 0; level--) {
			path[level] = node;
			node = node.parent;
		}

		double[] levelWeights = new double[numLevels];

		// Initialize level counts
		for (token = 0; token < seqLen; token++) {
			levelCounts[ docLevels[token] ]++;
		}

		for (token = 0; token < seqLen; token++) {
			type = fs.getIndexAtPosition(token);
	    
			levelCounts[ docLevels[token] ]--;
			node = path[ docLevels[token] ];
			node.typeCounts[type]--;
			node.totalTokens--;

			sum = 0.0;
			for (level=0; level < numLevels; level++) {
				double d = ((levelCounts[level]+1))*(level+1);
				if(Double.isInfinite(d) || Double.isNaN(d))
					throw new RuntimeException("Unclear level weight "+d);
				
				double e = getEta(type, level);
				double s = getEtaSum(type, level);
				double r = e/s;
				
				// guard against division by 0
				if(Double.isNaN(r))
					r = 0;
				
				levelWeights[level] = r;

				sum += levelWeights[level];
			}
			level = random.nextDiscrete(levelWeights, sum);

			docLevels[token] = level;
			levelCounts[ docLevels[token] ]++;
			node = path[ level ];
			node.typeCounts[type]++;
			node.totalTokens++;
		}
    }
    
    
    
    private double getEta(int type, int level)
    {
    	if(type >= dfs.length)
    		throw new RuntimeException("Type "+type+" is out of bounds "+dfs.length);

    	int bucket = (int)Math.floor(Utils.logb(dfs[type]+1, base));
    	
    	if(bucket < 0 || bucket >= numLevels)
    		throw new RuntimeException("Insane bucket "+bucket+" for "+numLevels+" levels");
    	    	
    			
    	if(bucket == level)
    	{
    		return eta;
    	}
    	else
    	{
    		return 0;
    	}
    }
    
    private double getEtaSum(int type, int level)
    {
    	return etaSum[level];
    }

    
	/**
	 *  Writes the current sampling state to the file specified in <code>stateFile</code>.
	 */
	public void printState() throws IOException, FileNotFoundException {
		printState(new PrintWriter(new BufferedWriter(new FileWriter(stateFile))));
	}

	/**
	 *  Write a text file describing the current sampling state. 
	 */
    public void printState(PrintWriter out) throws IOException {
		int doc = 0;

		Alphabet alphabet = instances.getDataAlphabet();

		for (Instance instance: instances) {
			FeatureSequence fs = (FeatureSequence) instance.getData();
			int seqLen = fs.getLength();
			int[] docLevels = levels[doc];
			NCRPNode node;
			int type, token, level;

			StringBuffer path = new StringBuffer();
			
			// Start with the leaf, and build a string describing the path for this doc
			node = documentLeaves[doc];
			for (level = numLevels - 1; level >= 0; level--) {
				path.append(node.nodeID + " ");
				node = node.parent;
			}

			for (token = 0; token < seqLen; token++) {
				type = fs.getIndexAtPosition(token);
				level = docLevels[token];
				
				// The "" just tells java we're not trying to add a string and an int
				out.println(path + "" + type + " " + alphabet.lookupObject(type) + " " + level + " ");
			}

			doc++;
		}
	}	    


    public void printNodes() {
		printNode(rootNode, 0);
    }

    public void printNode(NCRPNode node, int indent) {
    	if(node.totalTokens == 0)
    		return;
		StringBuffer out = new StringBuffer();
		for (int i=0; i<indent; i++) {
			out.append("  ");
		}

		
		out.append(node.getTopWords(numWordsToDisplay));
		System.out.println(out);
	
		for (NCRPNode child: node.children) {
			printNode(child, indent + 1);
		}
    }
    

    /** For use with empirical likelihood evaluation: 
     *   sample a path through the tree, then sample a multinomial over
     *   topics in that path, then return a weighted sum of words.
     */
    public double empiricalLikelihood(int numSamples, InstanceList testing)  {
		NCRPNode[] path = new NCRPNode[numLevels];
		NCRPNode node;
		double weight;
		path[0] = rootNode;

		FeatureSequence fs;
		int sample, level, type, token, doc, seqLen;

		Dirichlet dirichlet = new Dirichlet(numLevels, alpha);
		double[] levelWeights;
		double[] multinomial = new double[numTypes];

		double[][] likelihoods = new double[ testing.size() ][ numSamples ];

		for (sample = 0; sample < numSamples; sample++) {
			Arrays.fill(multinomial, 0.0);

			for (level = 1; level < numLevels; level++) {
				path[level] = path[level-1].selectExisting();
			}
	    
			levelWeights = dirichlet.nextDistribution();
	    
			for (type = 0; type < numTypes; type++) {
				for (level = 0; level < numLevels; level++) {
					node = path[level];
					multinomial[type] +=
						
//							levelWeights[level] *
						
//						(getEta(type, level) + node.typeCounts[type]) /
//						(getEtaSum(type, level) + node.totalTokens);
							(phi + node.typeCounts[type]) /
							(phiSum + node.totalTokens);
				}

			}

			for (type = 0; type < numTypes; type++) {
				multinomial[type] = Math.log(multinomial[type]);
			}

			for (doc=0; doc<testing.size(); doc++) {
                fs = (FeatureSequence) testing.get(doc).getData();
                seqLen = fs.getLength();
                
                for (token = 0; token < seqLen; token++) {
                    type = fs.getIndexAtPosition(token);
                    likelihoods[doc][sample] += multinomial[type];
                }
            }
		}
	
        double averageLogLikelihood = 0.0;
        double logNumSamples = Math.log(numSamples);
        for (doc=0; doc<testing.size(); doc++) {
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


    class NCRPNode {
		int customers;
		ArrayList<NCRPNode> children;
		NCRPNode parent;
		int level;

		int totalTokens;
		int[] typeCounts;

		public int nodeID;

		public NCRPNode(NCRPNode parent, int dimensions, int level) {
			customers = 0;
			this.parent = parent;
			children = new ArrayList<NCRPNode>();
			this.level = level;

			//System.out.println("new node at level " + level);
	    
			totalTokens = 0;
			typeCounts = new int[dimensions];

			nodeID = totalNodes;
			totalNodes++;
		}

		public NCRPNode(int dimensions) {
			this(null, dimensions, 0);
		}

		public NCRPNode addChild() {
			NCRPNode node = new NCRPNode(this, typeCounts.length, level + 1);
			children.add(node);
			return node;
		}

		public boolean isLeaf() {
			return level == numLevels - 1;
		}

		public NCRPNode getNewLeaf() {
			NCRPNode node = this;
			for (int l=level; l<numLevels - 1; l++) {
				node = node.addChild();
			}
			return node;
		}

		public void dropPath() {
			NCRPNode node = this;
			node.customers--;
			if (node.customers == 0) {
				node.parent.remove(node);
			}
			for (int l = 1; l < numLevels; l++) {
				node = node.parent;
				node.customers--;
				if (node.customers == 0) {
					node.parent.remove(node);
				}
			}
		}

		public void remove(NCRPNode node) {
			children.remove(node);
		}

		public void addPath() {
			NCRPNode node = this;
			node.customers++;
			for (int l = 1; l < numLevels; l++) {
				node = node.parent;
				node.customers++;
			}
		}

		public NCRPNode selectExisting() {
			double[] weights = new double[children.size()];
	    
			int i = 0;
			for (NCRPNode child: children) {
				weights[i] = (double) child.customers / (gamma + customers);
				i++;
			}

			int choice = random.nextDiscrete(weights);
			return children.get(choice);
		}

		public NCRPNode select() {
			double[] weights = new double[children.size() + 1];
	    
			weights[0] = gamma / (gamma + customers);

			int i = 1;
			for (NCRPNode child: children) {
				weights[i] = (double) child.customers / (gamma + customers);
				i++;
			}

			int choice = random.nextDiscrete(weights);
			if (choice == 0) {
				return(addChild());
			}
			else {
				return children.get(choice - 1);
			}
		}
	
		double [] vect = null;
		public double [] getTopicVector() {
			if(vect == null)
			{
				vect = new double [numTypes];
				double sum = 0;
				for (int type=0; type < numTypes; type++) {
					sum += typeCounts[type];
				}
			
				for (int type=0; type < numTypes; type++) {
					vect[type] = (typeCounts[type] + 0.05)/(sum + (0.05*normVect.length));
				}
			}
			return vect;
		}
		
		public double getExpectedRank()
		{
			double [] vect = getTopicVector();
			
			
			double exp = 0;
			for (int type=0; type < numTypes; type++) {
				exp += (vect[type] * dfs[type] );
			}
			return exp;
		}
		
		public String getTopWords(int numWords) {
			IDSorter[] sortedTypes = new IDSorter[numTypes];
	    
			for (int type=0; type < numTypes; type++) {
				sortedTypes[type] = new IDSorter(type, typeCounts[type]);
			}
			Arrays.sort(sortedTypes);
	    
			Alphabet alphabet = instances.getDataAlphabet();
			StringBuffer out = new StringBuffer();
			for (int i=0; i<Math.min(10, sortedTypes.length); i++) {
				if(typeCounts[sortedTypes[i].getID()] > 0)
				{
					out.append(alphabet.lookupObject(sortedTypes[i].getID()) +" ");
//					out.append(alphabet.lookupObject(sortedTypes[i].getID()) +"["+typeCounts[sortedTypes[i].getID()]+"]" + " ");
				}
			}
			return out.toString();
		}

    }
}

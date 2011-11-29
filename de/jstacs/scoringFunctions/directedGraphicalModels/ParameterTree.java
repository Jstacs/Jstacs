/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Jstacs is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.scoringFunctions.directedGraphicalModels;

import java.util.Arrays;
import java.util.LinkedList;

import de.jstacs.NonParsableException;
import de.jstacs.Storable;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sequence;
import de.jstacs.io.XMLParser;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.random.DirichletMRG;
import de.jstacs.utils.random.DirichletMRGParams;
import de.jtem.numericalMethods.calculus.specialFunctions.Gamma;

/**
 * Class for the tree that represents the context of a {@link Parameter} in a
 * {@link BayesianNetworkScoringFunction}.
 * 
 * @author Jan Grau
 */
public class ParameterTree implements Cloneable, Storable {

	private int pos;
	private int[] contextPoss;
	private TreeElement root;
	private AlphabetContainer alphabet;
	private int firstParent;
	private int[] firstChildren;

	/**
	 * Creates a new {@link ParameterTree} for the parameters at position
	 * <code>pos</code> using the parent positions in <code>contextPoss</code>.
	 * These are used to extract the correct alphabet out of
	 * <code>alphabet</code> for every context position. The first parent is the
	 * first parent of the random variable at <code>pos</code> as given by the
	 * topological ordering of the network structure of the enclosing
	 * {@link BayesianNetworkScoringFunction}. The first children are the
	 * children the random variable at <code>pos</code> is the first parent for.
	 * 
	 * @param pos    the position of the random variable of the parameters in the
	 *            tree
	 * @param contextPoss  the positions of the context
	 * @param alphabet  the alphabet of the enclosing
	 *            {@link BayesianNetworkScoringFunction}
	 * @param firstParent   the first parent of this random variable, or <code>-1</code>
	 *            if the random variable has no parent
	 * @param firstChildren  the first children of this random variable
	 */
	public ParameterTree(int pos, int[] contextPoss,
			AlphabetContainer alphabet, int firstParent, int[] firstChildren) {
		this.pos = pos;
		this.contextPoss = contextPoss;
		this.alphabet = alphabet;
		this.firstParent = firstParent;
		this.firstChildren = firstChildren;
		this.root = new TreeElement(0, alphabet);
	}

	/**
	 * Recreates a {@link ParameterTree} from its XML representation as
	 * returned by {@link #toXML()}. {@link ParameterTree} does not
	 * implement the {@link Storable} interface to recycle the
	 * {@link AlphabetContainer} of the enclosing
	 * {@link BayesianNetworkScoringFunction}, but besides the different
	 * constructor works like any implementation of {@link Storable}.
	 * 
	 * @param source  the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException if the XML code could not be parsed
	 */
	public ParameterTree(StringBuffer source)
			throws NonParsableException {
		try {
			source = XMLParser.extractForTag(source, "parameterTree");
			pos = XMLParser.extractObjectForTags(source, "pos", int.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
			contextPoss = XMLParser.extractObjectForTags(source, "contextPoss", int[].class );
			root = new TreeElement(XMLParser.extractForTag(source, "root"));
			this.alphabet = null;
			this.firstParent = XMLParser.extractObjectForTags(source, "firstParent", int.class );
			this.firstChildren = XMLParser.extractObjectForTags(source, "firstChildren", int[].class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
		} catch (NonParsableException e) {
			e.printStackTrace();
			throw e;
		}

	}

	void setAlphabet(AlphabetContainer alphabet){
		this.alphabet = alphabet;
	}
	
	/* (non-Javadoc)
	 * @see java.lang.Object#clone()
	 */
	@Override
	public ParameterTree clone() throws CloneNotSupportedException {
		ParameterTree clone = (ParameterTree) super.clone();
		clone.contextPoss = contextPoss.clone();
		clone.cloneRoot();
		clone.firstChildren = firstChildren.clone();
		return clone;
	}

	private void cloneRoot() throws CloneNotSupportedException {
		TreeElement temp = this.root;
		this.root = new TreeElement(root.contNum, alphabet);
		this.root.cloneRest(temp);
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		StringBuffer all = new StringBuffer();
		all.append("Probabilities at position " + pos + ":\n");
		root.appendToBuffer(all, "");
		return all.toString();
	}

	/**
	 * Computes the probabilities for a PWM, i.e. the parameters in the tree
	 * have an empty context, and inserts them into <code>probs</code>. Used by
	 * {@link BayesianNetworkScoringFunction#getPWM()}.
	 * 
	 * @param probs  the array to store the probabilities for a PWM
	 * 
	 * @throws Exception if the tree structure does have a non-empty context
	 */
	public void insertProbs(double[] probs) throws Exception {
		root.insertProbs(probs);
	}

	/**
	 * Extracts the {@link Parameter}s from the leaves of this tree in
	 * left-to-right order (as specified by the order of the alphabet) and
	 * returns them as a {@link LinkedList}.
	 * 
	 * @return the {@link Parameter}s from the leaves in linear order
	 */
	public LinkedList<Parameter> linearizeParameters() {
		return root.linearizeParameters(new LinkedList<Parameter>());
	}

	/**
	 * Indicates if the random variable of this {@link ParameterTree} is a
	 * leaf, i.e. it has no children in the network structure of the enclosing
	 * {@link BayesianNetworkScoringFunction}.
	 * 
	 * @return <code>true</code> if this tree is a leaf, <code>false</code> otherwise
	 */
	public boolean isLeaf() {
		return firstChildren.length == 0;
	}

	/**
	 * Returns the number of parents for the random variable of this
	 * {@link ParameterTree} in the network structure of the enclosing
	 * {@link BayesianNetworkScoringFunction}. This corresponds to the length of
	 * the context or the depth of the tree.
	 * 
	 * @return the number of parents
	 */
	public int getNumberOfParents() {
		return contextPoss.length;
	}

	/**
	 * Prints the structure of this tree.
	 */
	public void print() {
		System.out.println("tree " + pos + ": ");
		root.print();
	}

	/**
	 * Returns the {@link Parameter} that is responsible for the suffix of
	 * sequence <code>seq</code> starting at position <code>start</code>.
	 * 
	 * @param seq the {@link Sequence}
	 * @param start  the first position in the suffix
	 * 
	 * @return the {@link Parameter} that is responsible for the suffix
	 */
	public Parameter getParameterFor(Sequence seq, int start) {
		return root.getParameterFor(seq, start);
	}

	/**
	 * Sets the instance of the {@link Parameter} for symbol <code>symbol</code> and
	 * context <code>context</code> to {@link Parameter} <code>par</code>.
	 * 
	 * @param symbol   the symbol
	 * @param context    the context
	 * @param par  the new {@link Parameter} instance
	 */
	public void setParameterFor(int symbol, int[][] context, Parameter par) {
		root.setParameterFor(0, symbol, context, par);
	}

	/**
	 * Resets all pre-computed normalization constants.
	 */
	public void invalidateNormalizers() {
		root.invalidateNormalizers();
	}

	/**
	 * Computes the forward-part of the normalization constant starting from
	 * this {@link ParameterTree}. This is only possible for roots, i.e.
	 * {@link ParameterTree}s that do not have parents in the network structure
	 * of the enclosing {@link BayesianNetworkScoringFunction}.
	 * 
	 * @param trees  the array of all trees as from the enclosing
	 *            {@link BayesianNetworkScoringFunction}
	 *            
	 * @return the forward-part of the normalization constant
	 * 
	 * @throws RuntimeException if this {@link ParameterTree} is not a root
	 */
	public double forward(ParameterTree[] trees) throws RuntimeException {
		if (this.getNumberOfParents() > 0) {
			throw new RuntimeException("Forward can only be started at roots.");
		} else {
			// System.out.println();
			return this.getLogZ(new int[0][2], trees);
		}
	}

	private double getLogZ(int[][] context, ParameterTree[] trees)
			throws RuntimeException {
		return root.getLogZ(context, new int[this.getNumberOfParents() + 1][2],
				trees, 0);
	}

	private double getLogT(int[][] context, ParameterTree[] trees, int[][] order)
			throws RuntimeException {
		return root.getLogT(
						context,
						firstParent > -1 ? new int[trees[firstParent].contextPoss.length + 1][2]
								: new int[0][2], trees, order, 0);
	}

	/**
	 * Computes the backward-part of the normalization constant
	 * starting from this {@link ParameterTree}. This is only possible for
	 * leaves, i.e. {@link ParameterTree}s that do not have children in the
	 * network structure of the enclosing {@link BayesianNetworkScoringFunction}.
	 * 
	 * @param trees the array of all trees as from the enclosing
	 *            {@link BayesianNetworkScoringFunction}
	 * @param order  the topological ordering as returned by
	 *            {@link de.jstacs.algorithms.graphs.TopSort#getTopologicalOrder(int[][])}
	 *            
	 * @throws RuntimeException if this {@link ParameterTree} is not a leaf
	 */
	public void backward(ParameterTree[] trees, int[][] order)
			throws RuntimeException {
		if (!this.isLeaf()) {
			throw new RuntimeException(
					"Backward can only be started at leaves.");
		} else {
			root.startBackward(new int[this.getNumberOfParents() + 1][2],
					trees, order, 0);
		}
	}

	/**
	 * Adds <code>count</code> to the parameter as returned by
	 * {@link #getParameterFor(Sequence, int)}.
	 * 
	 * @param seq  the sequence
	 * @param start  the first position of the suffix of <code>seq</code>
	 * @param count   the added count
	 */
	public void addCount(Sequence seq, int start, double count) {
		this.getParameterFor(seq, start).addCount(count);
	}

	/**
	 * Starts the normalization of the plug-in parameters to the logarithm of the
	 * MAP-estimates.
	 */
	public void normalizePlugInParameters() {
		root.normalizePlugInParameters();
	}

	/**
	 * Normalizes the parameter values to the corresponding log-probabilities.
	 * After this step,
	 * {@link BayesianNetworkScoringFunction#getLogNormalizationConstant()} should
	 * return 0.
	 */
	public void normalizeParameters() {
		root.normalizeParameters();
	}

	/**
	 * Divides each of the normalized parameters on a simplex by the last
	 * {@link Parameter}, which is defined not to be free. In the log-space this amounts
	 * to subtracting the value of the last {@link Parameter} from all {@link Parameter}s on the
	 * simplex.
	 */
	public void divideByUnfree() {
		root.divideByUnfree();
	}

	/**
	 * Works as defined in {@link Storable}. 
	 * Returns an XML representation of this {@link ParameterTree}.
	 * 
	 * @return the XML representation of this {@link ParameterTree}
	 * 
	 * @see Storable#toXML()
	 */
	public StringBuffer toXML() {
		StringBuffer source = new StringBuffer();
		XMLParser.appendObjectWithTags(source, pos, "pos");
		XMLParser.appendObjectWithTags(source, contextPoss, "contextPoss");
		XMLParser.appendObjectWithTags(source, root, "root");
		XMLParser.appendObjectWithTags(source, firstParent, "firstParent");
		XMLParser
				.appendObjectWithTags(source, firstChildren, "firstChildren");
		XMLParser.addTags(source, "parameterTree");
		return source;

	}

	/*
	 * Creates a {@link String} array from the XML representation of all
	 * {@link ParameterTree}s in <code>trees</code> as returned by the
	 * corresponding {@link #toXML()} methods. 
	 * 
	 * @param trees   the array of {@link ParameterTree}s, e.g. as returned by {@link #fromStringArray(String[], AlphabetContainer)}
	 * 
	 * @return the XML representations of the {@link ParameterTree}s as {@link String} array
	 * 
	 * @see ParameterTree#fromStringArray(String[], AlphabetContainer)
	 */
	/*public static String[] toStringArray(ParameterTree[] trees) {
		String[] strs = new String[trees.length];
		for (int i = 0; i < trees.length; i++) {
			strs[i] = trees[i].toXML().toString();
		}
		return strs;
	}*/

	/*
	 * Recreates an array of {@link ParameterTree}s from their
	 * XML representation as given in <code>strs</code>.
	 * 
	 * @param strs the XML representations as {@link String} array, e.g. as returned by
	 *            {@link #toStringArray(ParameterTree[])}
	 * @param alphabet  the alphabet of the enclosing
	 *            {@link BayesianNetworkScoringFunction}
	 *            
	 * @return the array of {@link ParameterTree}s
	 * 
	 * @throws NonParsableException if one of the XML representations could not be
	 *             parsed
	 *  
	 * @see ParameterTree#toStringArray(ParameterTree[]) 
	 */
	/*public static ParameterTree[] fromStringArray(String[] strs,
			AlphabetContainer alphabet) throws NonParsableException {
		ParameterTree[] trees = new ParameterTree[strs.length];
		for (int i = 0; i < trees.length; i++) {
			trees[i] = new ParameterTree(new StringBuffer(strs[i]), alphabet);
		}
		return trees;
	}*/

	/**
	 * Returns the first parent of the random variable of this
	 * {@link ParameterTree} in the topological ordering of the network
	 * structure of the enclosing {@link BayesianNetworkScoringFunction}.
	 * 
	 * @return the first parent
	 */
	public int getFirstParent() {
		return firstParent;
	}

	/**
	 * Draws KL-divergences between the distribution given by <code>contrast</code> and 
	 * <code>endIdx-startIdx</code> distributions drawn from a Dirichlet density centered around <code>contrast</code>, i.e. the hyper-parameters
	 * of the Dirichlet density are the probabilities of <code>contrast</code> weighted by <code>samples</code>. The drawn KL-divergences are stored
	 * into <code>kls</code> between <code>startIndex</code> and <code>endIndex</code> (exclusive).
	 * 
	 * @param weight a weight on the KL-divergences
	 * @param kls the array of KL-divergences which is filled
	 * @param startIdx the first index
	 * @param endIdx the index after the last index
	 * @param contrast the distribution to check against
	 * @param samples number of sequences + equivalent sample size
	 */
	public void drawKLDivergences(double weight, double[] kls, int startIdx, int endIdx,
			double[][][] contrast, double samples) {
		root.drawKLDivergences(weight, kls, startIdx, endIdx, contrast, samples, 0, 0);
	}

	/**
	 * Returns the KL-divergence of the distribution of this {@link ParameterTree} and the distribution given by
	 * <code>ds</code>.
	 * 
	 * @param ds the distribution
	 * @return the KL-divergence
	 */
	public double getKLDivergence(double[][][] ds) {
		return root.getWeightedKLDivergence(ds, 0, 0);
	}
	
	/**
	 * Returns the KL-divergence of the distribution of this {@link ParameterTree} and a number of distribution given by
	 * <code>ds</code> and weighted by <code>weight</code>
	 * 
	 * @param distribution the distribution
	 * @param weight the weights on the distributions
	 * @return the KL-divergence
	 */
	public double getKLDivergence(double[] weight, double[][][][] distribution ) {
		return root.getWeightedKLDivergence( weight, distribution, 0, 0 );
	}
	
	/**
	 * Draws KL-divergences between the distributions given by <code>contrast[i]</code> each weighted by <code>weights[i]</code> 
	 * <code>kls.length</code> distributions drawn from a Dirichlet density centered around <code>contrast</code>, i.e. the hyper-parameters
	 * of the Dirichlet density are the probabilities of <code>contrast</code> weighted by <code>samples</code>. The drawn KL-divergences are added
	 * to the entries of <code>kls</code>.
	 * 
	 * @param kls the array of KL-divergences which is filled
	 * @param weights the weights on the distributions in contrast
	 * @param contrast the distribution to check against
	 * @param samples number of sequences + equivalent sample size
	 */
	public void drawKLDivergences(double[] kls, double[] weights, double[][][][] contrast, double samples) {
		root.drawKLDivergences( kls, weights, contrast, samples, 0, 0 );
	}

	/**
	 * Fills all parameters with the probabilities given in
	 * <code>distribution</code>.
	 * 
	 * @param distribution the distribution
	 * @param weight the weights on the distributions
	 */
	public void fill(double[] weight, double[][][][] distribution) {
		
		root.setNewParameters( weight, distribution, 0, 0 );
	}

	/**
	 * Copies the values of the parameters from another {@link ParameterTree}.
	 * 
	 * @param parameterTree   the template
	 */
	public void copy(ParameterTree parameterTree) {
		root.copy(parameterTree.root);
	}

	/**
	 * Initializes the parameters of this {@link ParameterTree} randomly.
	 * 
	 * @param ess    the equivalent sample size
	 */
	public void initializeRandomly(double ess) {
		root.initializeRandomly(ess);
	}

	/**
	 * Computes the Gamma-normalization for the prior.
	 * 
	 * @return the Gamma-normalization
	 */
	public Double computeGammaNorm() {
		return root.computeGammaNorm();
	}
	
	/**
	 * Returns the probability of {@link Sequence} <code>sequence</code> in this {@link ParameterTree}.
	 * @param sequence the sequence
	 * @return the probability
	 */
	public double getProbFor( Sequence sequence ) {
		return root.getProbFor(sequence, 0);
	}

	/**
	 * Class for the nodes of a {@link ParameterTree}
	 * @author Jan Grau
	 *
	 */
	public class TreeElement implements Storable, Cloneable {

		private int contextPos;
		private TreeElement[] children;
		private Parameter[] pars;
		private Double fullNormalizer;
		private Double[] symT;
		private int contNum;

		private TreeElement(int contNum, AlphabetContainer alphabet) {
			this.contNum = contNum;
			if (contNum < contextPoss.length) {
				this.contextPos = contextPoss[contNum];
				children = new TreeElement[(int) alphabet
						.getAlphabetLengthAt(this.contextPos)];
				for (int i = 0; i < alphabet
						.getAlphabetLengthAt(this.contextPos); i++) {
					children[i] = new TreeElement(contNum + 1, alphabet);
				}
			} else {
				this.contextPos = -1;
				this.pars = new Parameter[(int) alphabet
						.getAlphabetLengthAt(pos)];

				this.fullNormalizer = null;
				this.symT = new Double[pars.length];
			}
		}

		private void appendToBuffer(StringBuffer all, String after) {
			if (children != null) {
				for (int i = 0; i < children.length; i++) {
					children[i].appendToBuffer(all, after + "X_" + contextPos
							+ " = " + alphabet.getSymbol(contextPos, i) + ", ");
				}
			} else {
				double[] norms = new double[pars.length];;
				for (int i = 0; i < pars.length; i++) {
					norms[i] = pars[i].getValue() + pars[i].getLogZ();
				}
				double logNorm = Normalisation.getLogSum( norms );
				for (int i = 0; i < pars.length; i++) {
					double tempTheta = Math.exp(pars[i].getValue() + pars[i].getLogZ()
							- logNorm);
					all.append("P(X_" + pos + " = "
							+ alphabet.getSymbol(pos, i) + " | " + after
							+ "c)=" + tempTheta);
					if (i < pars.length - 1) {
						all.append("\t");
					}
				}
				all.append("\n");
			}
		}

		private void normalizeParameters() {
			if (children != null) {
				for (int i = 0; i < children.length; i++) {
					children[i].normalizeParameters();
				}
			} else {
				double[] norms = new double[pars.length];;
				for (int i = 0; i < pars.length; i++) {
					norms[i] = pars[i].getValue() + pars[i].getLogZ();
				}
				double logNorm = Normalisation.getLogSum( norms );
				for (int i = 0; i < pars.length; i++) {
					double tempTheta = pars[i].getValue() + pars[i].getLogZ()
							- logNorm;
					pars[i].setValue(tempTheta);
				}
				if (!pars[pars.length - 1].isFree()) {
					for (int i = 0; i < pars.length; i++) {
						pars[i].setValue(pars[i].getValue()
								- pars[pars.length - 1].getValue());
					}
				}
			}
		}

		private void insertProbs(double[] probs) throws Exception {
			if (children != null) {
				throw new Exception("Not implemented");
			} else {
				for (int i = 0; i < pars.length; i++) {
					probs[i] = pars[i].getValue() + pars[i].getLogZ();
				}
				Normalisation.logSumNormalisation(probs);

			}
		}

		private void startBackward(int[][] newContext, ParameterTree[] trees,
				int[][] order, int depth) throws RuntimeException {
			if (children != null) {
				newContext[depth][0] = contextPos;
				for (int i = 0; i < children.length; i++) {
					newContext[depth][1] = i;
					children[i].startBackward(newContext, trees, order,
							depth + 1);
				}
			} else {
				newContext[depth][0] = pos;
				for (int i = 0; i < pars.length; i++) {
					newContext[depth][1] = pars[i].symbol;
					trees[pos].getLogT(newContext, trees, order);
				}
			}
		}

		private double getLogT(int[][] context, int[][] newContext,
				ParameterTree[] trees, int[][] order, int depth)
				throws RuntimeException {
			// we are not at a leaf of the ParameterTree
			if (children != null) {
				// search for my context's position
				for (int i = 0; i < context.length; i++) {
					if (context[i][0] == this.contextPos) {
						// copy to newContext, we might need it for the
						// predecessors, if the have the same
						newContext[depth][0] = context[i][0];
						newContext[depth][1] = context[i][1];
						// go to child for found context
						return children[context[i][1]].getLogT(context,
								newContext, trees, order, depth + 1);
					}
				}
				// should not happen
				throw new RuntimeException(
						"Correct context not found for depth " + depth
								+ " at position " + pos + ".");
			} else {
				// search for the context that corresponds to my position
				for (int i = 0; i < context.length; i++) {
					if (context[i][0] == pos) {
						// search for parameter that corresponds to realization
						// of found context
						for (int j = 0; j < pars.length; j++) {
							if (pars[j].symbol == context[i][1]) {
								// if already computed for this parameter, done
								if (symT[j] != null) {
									return symT[j];
								} else {
									int fp = firstParent;
									// if we are at a root
									if (fp == -1) {
										pars[j].setLogT(0d);
										return pars[j].getValue();
										// if we are at a node with a context
										// that fully defines its first parent,
										// i.e. the parent that fully defines
										// the configuration of this node
									} else if (trees[fp].contextPoss.length < contextPoss.length) {
										// recursion
										double temp = trees[fp].getLogT(
												newContext, trees, order);
										int[] fcoffp = trees[fp].firstChildren;
										// for the other branches under fp that
										// are independent of this node's
										// branch, we need the Zs, which have
										// already been computed
										for (int k = 0; k < fcoffp.length; k++) {
											// but not for this branch
											if (fcoffp[k] != pos) {
												// System.out.println(fp+" has fc "+fcoffp[k]);
												temp += trees[fcoffp[k]].getLogZ(
														newContext, trees);
											}
										}
										// set the T for the partial
										// normalization
										pars[j].setLogT(temp);
										// compute the T for recursion
										symT[j] = pars[j].getValue() + temp;
										return symT[j];
										// we are at a node with a parent that
										// has more and different parents than
										// this node
										// (at most one, this is the parent over
										// whose configs we will sum)
									} else {

										// find the parent of the first parent
										// that has the lowest number in the
										// topological order,
										// this must be the parent on which we
										// do not depend
										int[] cp = trees[fp].contextPoss;
										int lowestOrder = Integer.MAX_VALUE;
										int lowestOrderIndex = -1;
										for (int k = 0; k < cp.length; k++) {
											if (order[cp[k]][1] < lowestOrder) {
												lowestOrder = order[cp[k]][1];
												lowestOrderIndex = cp[k];
											}
										}

										// sum over configurations of the parent
										// lowestOrderIndex
										newContext[depth][0] = lowestOrderIndex;
										int al = (int) alphabet.getAlphabetLengthAt(lowestOrderIndex);
										double[] temp = new double[al];
										for (byte a = 0; a < al; a++) {
											newContext[depth][1] = a;
											// recursion
											temp[a] = trees[fp].getLogT(
													newContext, trees, order);
											// for the other branches under fp
											// that are independent of this
											// node's
											// branch, we need the Zs, which
											// have already been computed
											int[] fcoffp = trees[fp].firstChildren;
											for (int k = 0; k < fcoffp.length; k++) {
												// but nor for this branch
												if (fcoffp[k] != pos) {
													temp[a] += trees[fcoffp[k]]
															.getLogZ(newContext,
																	trees);
												}
											}
										}
										double temp2 = Normalisation.getLogSum( temp );
										// set T for partial normalization
										pars[j].setLogT(temp2);
										// compute the T for recursion
										symT[j] = pars[j].getValue() + temp2;
										return symT[j];
									}
								}
							}
						}
					}
				}
				throw new RuntimeException(
						"Parameter value not defined in context.");
			}
		}

		private double getLogZ(int[][] context, int[][] newContext,
				ParameterTree[] trees, int depth) throws RuntimeException {
			// we are not at a leaf of the ParameterTree
			if (children != null) {
				// search for my context's position
				for (int i = 0; i < context.length; i++) {
					// found
					if (context[i][0] == this.contextPos) {
						// copy to newContext, we need it for the descendants
						newContext[depth][0] = context[i][0];
						newContext[depth][1] = context[i][1];
						// go to children for found contet
						return children[context[i][1]].getLogZ(context,
								newContext, trees, depth + 1);
					}
				}
				// should not happen
				throw new RuntimeException(
						"Correct context could not be found at position " + pos
								+ " and depth " + depth);
			} else if (fullNormalizer != null) {
				// we already precomputed the normalization constant, done
				return fullNormalizer;
			} else {
				double[] vals = new double[pars.length];
				// to compute the normalization constant we must sum over all
				// children
				for (int i = 0; i < pars.length; i++) {
					int[] fc = firstChildren;
					if (fc == null) {
						throw new RuntimeException(
								"First children of parameter "
										+ pars[i].getIndex() + " not defined.");
					}
					// this parameter's value will be the context of its
					// descendants
					newContext[depth][0] = pars[i].getPosition();
					newContext[depth][1] = pars[i].symbol;
					double temp = 0;
					// compute product over all first children, i.e. children
					// that are the roots
					// of independent subtrees and whose context is fully
					// defined by the current parameter
					for (int j = 0; j < fc.length; j++) {
						temp += trees[fc[j]].getLogZ(newContext, trees);
					}
					// set the Z-part of the local normalization (for partial
					// normalization constant)
					pars[i].setLogZ(temp);
					// for fullNormalizer multiplied by current value
					vals[i] = pars[i].getValue() + temp;
				}
				fullNormalizer = Normalisation.getLogSum(vals);
				return fullNormalizer;
			}
		}

		private void invalidateNormalizers() {
			if (children != null) {
				for (int i = 0; i < children.length; i++) {
					children[i].invalidateNormalizers();
				}
			} else {
				for (int i = 0; i < pars.length; i++) {
					pars[i].invalidateNormalizers();
					symT[i] = null;
				}
			}
			fullNormalizer = null;
		}

		private void cloneRest(TreeElement original)
				throws CloneNotSupportedException {
			contextPos = original.contextPos;
			if (children != null) {
				this.children = new TreeElement[children.length];
				for (int i = 0; i < children.length; i++) {
					this.children[i] = new TreeElement(
							original.children[i].contNum, alphabet);
					children[i].cloneRest(original.children[i]);
				}
			} else {
				children = null;
			}
			if (pars != null) {
				pars = new Parameter[pars.length];
				for (int i = 0; i < pars.length; i++) {
					pars[i] = original.pars[i].clone();
				}
				fullNormalizer = null;
				symT = new Double[pars.length];
			}
		}

		private void divideByUnfree() {
			if (pars != null) {
				double div = pars[pars.length - 1].getValue();
				for (int i = 0; i < pars.length; i++) {
					if (!Double.isNaN(pars[i].getValue() - div)
							&& !Double.isInfinite(pars[i].getValue() - div)) {
						pars[i].setValue(pars[i].getValue() - div);
					} else {
						pars[i].setValue(0d);
					}
				}
			} else {
				for (int i = 0; i < children.length; i++) {
					children[i].divideByUnfree();
				}
			}
		}

		private LinkedList<Parameter> linearizeParameters(
				LinkedList<Parameter> list) {
			if (children != null) {
				for (int i = 0; i < children.length; i++) {
					children[i].linearizeParameters(list);
				}
			} else {
				for (int i = 0; i < pars.length; i++) {
					list.add(pars[i]);
				}
			}
			return list;
		}

		/**
		 * Constructor for the {@link Storable} interface.
		 * @param representation the XML-representation of the {@link TreeElement}
		 * @throws NonParsableException if the XML could not be parsed
		 */
		public TreeElement(StringBuffer representation)
				throws NonParsableException {
			representation = XMLParser.extractForTag(representation,
					"treeElement");
			contNum = XMLParser.extractObjectForTags(representation, "contNum", int.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
			contextPos = XMLParser.extractObjectForTags(representation, "contextPos", int.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
			
			//XXX not possible since inner class
			children = XMLParser.extractObjectAndAttributesForTags( representation, "children", null, null, TreeElement[].class, ParameterTree.class, ParameterTree.this );
			/*
			StringBuffer temp = XMLParser.extractForTag(representation,"children");
			if (temp.toString().equals("null")) {
				children = null;
			} else {
				XMLParser.addTags(temp, "children");
				String[] childRep = XMLParser.extractObjectForTags(temp, "children", String[].class );
				children = new TreeElement[childRep.length];
				for (int i = 0; i < childRep.length; i++) {
					children[i] = new TreeElement(new StringBuffer(childRep[i]));
				}
			}
			*/
			
			pars = XMLParser.extractObjectForTags(representation, "pars", Parameter[].class );
			/*
			temp = XMLParser.extractForTag(representation, "pars");
			if (temp.toString().equals("null")) {
				pars = null;
			} else {
				XMLParser.addTags(temp, "pars");
				pars = XMLParser.extractObjectForTags(temp, "pars", Parameter[].class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
			}
			*/
			if (pars != null) {
				symT = new Double[pars.length];
				fullNormalizer = null;
			}
		}

		private void setParameterFor(int depth, int symbol, int[][] context,
				Parameter par) {
			if (children != null) {
				// System.out.println("going to child for "+context[depth][1]+" in tree "+pos);
				for (int i = 1; i < context[depth].length; i++) {
					children[context[depth][i]].setParameterFor(depth + 1,
							symbol, context, par);
				}
			} else {
				// System.out.println("setting parameter for "+symbol);
				pars[symbol] = par;
			}
		}

		private void print() {
			System.out.println(contextPos);
			if (children != null) {
				for (int i = 0; i < children.length; i++) {
					System.out.println("child " + i + ":");
					children[i].print();
				}
			} else {
				for (int i = 0; i < pars.length; i++) {
					pars[i].print();
				}
			}
		}

		private void normalizePlugInParameters() {
			if (children != null) {
				for (int i = 0; i < children.length; i++) {
					children[i].normalizePlugInParameters();
				}
			} else {
				double sum = 0;
				for (int i = 0; i < pars.length; i++) {
					sum += pars[i].getCounts();
				}
				if (sum > 0) {
					for (int i = 0; i < pars.length; i++) {
						pars[i].setValue(Math.log(pars[i].getCounts() / sum));
					}
				} else {
					for (int i = 0; i < pars.length; i++) {
						pars[i].setValue(-Math.log(pars.length));
					}
				}
			}
		}

		private Parameter getParameterFor(Sequence seq, int start) {
			if (children != null) {
				return children[seq.discreteVal(contextPos + start)]
						.getParameterFor(seq, start);
			} else {
				return pars[seq.discreteVal(pos + start)];
			}
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see de.jstacs.Storable#toXML()
		 */
		public StringBuffer toXML() {
			StringBuffer source = new StringBuffer();
			XMLParser.appendObjectWithTags(source, contNum, "contNum");
			XMLParser.appendObjectWithTags(source, contextPos, "contextPos");
			
			XMLParser.appendObjectWithTags(source, children, "children");
			XMLParser.appendObjectWithTags(source, pars, "pars");
			/*
			if (children != null) {
				XMLParser.appendObjectWithTags(source, children,
						"children");
			} else {
				XMLParser.appendObjectWithTags(source, "null", "children");
			}
			if (pars != null) {
				XMLParser.appendObjectWithTags(source, pars, "pars");
			} else {
				XMLParser.appendObjectWithTags(source, "null", "pars");
			}
			*/
			XMLParser.addTags(source, "treeElement");
			return source;
		}

		private void drawKLDivergences(double weight, double[] kls, int startIdx, int endIdx,
				double[][][] ds, double samples, int context, int depth) {
			if (children != null && depth < ds.length - 1) {
				for (int i = 0; i < children.length; i++) {
					children[i].drawKLDivergences(weight, kls, startIdx, endIdx, ds,
							samples, context + i * ds[depth].length, depth + 1);
				}
			} else if (children != null) {
				for (int i = 0; i < children.length; i++) {
					children[i].drawKLDivergences(weight, kls, startIdx, endIdx, ds,
							samples, context, depth);
				}
			} else {
				double[] dist = ds[depth][context];
				double[] weightedEss = new double[this.pars.length];
				double contextProb = getContextProbability();
				for (int i = 0; i < weightedEss.length; i++) {
					weightedEss[i] = samples * contextProb * dist[i];
				}

				DirichletMRGParams p = new DirichletMRGParams(weightedEss);
				double[] temp = new double[weightedEss.length];
				for (int i = startIdx; i < endIdx; i++) {
					DirichletMRG.DEFAULT_INSTANCE.generate(temp, 0,
							temp.length, p);
					for (int j = 0; j < temp.length; j++) {
						if (temp[j] > 0.0) {
							kls[i] += weight * contextProb * temp[j]
									* Math.log(temp[j] / dist[j]);
						}
					}
				}

			}
		}
		
		private void setNewParameters(double[] weight, double[][][][] distribution, int context,
				int depth) {

			int i = 0, a = (int) alphabet.getAlphabetLengthAt( contextPos );
			if( children != null ) {
				a = (int) Math.pow( a, depth );
				for( ; i < children.length; i++ ) {
					children[i].setNewParameters(weight, distribution, context + i	* a, depth + 1);
				}
			} else {

				//fill in marginal distribution
				fill(getMarginal( weight, distribution, context, depth ));
			}
		}
		

		private double getWeightedKLDivergence(double[][][] ds, int context,
				int depth) {
			double kl = 0;
			if (children != null && depth < ds.length - 1) {

				for (int i = 0; i < children.length; i++) {
					kl += children[i].getWeightedKLDivergence(ds, context + i
							* ds[depth].length, depth + 1);
				}
			} else if (children != null) {
				for (int i = 0; i < children.length; i++) {
					kl += children[i].getWeightedKLDivergence(ds, context,
							depth);
				}
			} else {
				double val = 0;
				double[] norms = new double[pars.length];;
				for (int i = 0; i < pars.length; i++) {
					norms[i] = pars[i].getValue() + pars[i].getLogZ();
				}
				double logNorm = Normalisation.getLogSum( norms );
				// System.out.print("comparing "+Arrays.toString(
				// ds[depth][context] )+" to [");
				for (int i = 0; i < pars.length; i++) {
					double temp = Math.exp(pars[i].getValue() + pars[i].getLogZ() - logNorm);
					// System.out.print(temp+", ");
					if (temp > 0) {
						val += temp * Math.log(temp / ds[depth][context][i]);
					}
				}
				// System.out.println("] => "+val);

				double weight = getContextProbability();
				kl = val * weight;
			}
			return kl;
		}
		
		private double getWeightedKLDivergence(double[] weight, double[][][][] distribution, int context,
				int depth) {
			double kl = 0;
			int i = 0, a = (int) alphabet.getAlphabetLengthAt( contextPos );
			if( children != null ) {
				a = (int) Math.pow( a, depth );
				for( ; i < children.length; i++ ) {
					kl += children[i].getWeightedKLDivergence(weight, distribution, context + i	* a, depth + 1);
				}
			} else {
				double[] norms = new double[pars.length];
				double temp;
				//compute norm
				for( ; i < pars.length; i++ ) {
					norms[i] += pars[i].getValue() + pars[i].getLogZ();
				}
				double logNorm = Normalisation.getLogSum( norms );
				//compute marginal distribution
				double[] temp2 = getMarginal( weight, distribution, context, depth );
				//compute KL between motif and marginal distribution
				for( i = 0; i < pars.length; i++ ) {
					temp = Math.exp(pars[i].getValue() + pars[i].getLogZ() - logNorm);
					if (temp > 0) {
						kl += temp * Math.log(temp / temp2[i]);
					}
				}
				kl *= getContextProbability();
			}
			return kl;
		}
		
		private double[] getMarginal(double[] weight, double[][][][] distribution, int context,
				int depth){
			int a = (int) alphabet.getAlphabetLengthAt( contextPos );
			//compute marginal distribution
			double[] marginal = new double[pars.length];
			for( int c,d,j,i = 0; i < weight.length; i++ ) {
				if( depth < distribution[i].length ) {
					d = depth;
					c = context;
				} else {
					d = distribution[i].length-1;
					c = context % (int)Math.pow( a, d );
				}
				for( j = 0; j < pars.length; j++ ){
					marginal[j] += weight[i] * distribution[i][d][c][j];
				}
			}
			return marginal;
		}
		
		private void drawKLDivergences( double[] kls, double[] weight, double[][][][] distribution, double samples, int context, int depth ) {
			int i = 0, j, c, d,  a = (int) alphabet.getAlphabetLengthAt( contextPos );
			if (children != null) {
				a = (int) Math.pow( a, depth );
				for( ; i < children.length; i++ ) {
					children[i].drawKLDivergences( kls, weight, distribution, samples, context + i * a, depth + 1);
				}
			} else {
				double[] weightedEss = new double[pars.length], marginal = new double[pars.length];
				double contextProb = getContextProbability();
				DirichletMRGParams[] p = new DirichletMRGParams[weight.length];
				//compute marginal and create hyper-parameters
				for( i = 0; i < weight.length; i++ ){
					if( depth < distribution[i].length ) {
						d = depth;
						c = context;
					} else {
						d = distribution[i].length-1;
						c = context % (int)Math.pow( a, d );
					}
					for( j = 0; j < pars.length; j++ ){
						marginal[j] += weight[i] * distribution[i][d][c][j];
						weightedEss[j] = samples * contextProb * weight[i] * distribution[i][d][c][j];
					}
					p[i] = new DirichletMRGParams(weightedEss);
				}

				double[] part = new double[pars.length], marginalDrawn = new double[pars.length];
				for( i = 0; i < kls.length; i++) {
					// draw and compute marginal
					Arrays.fill(  marginalDrawn, 0 );
					for( j = 0; j < p.length; j++ ) {
						DirichletMRG.DEFAULT_INSTANCE.generate( part, 0, part.length, p[j] );
						for( a = 0; a < pars.length; a++ ) {
							marginalDrawn[a] += weight[j] * part[a];
						}
					}
					// compute the i-th KL
					for( j = 0; j < part.length; j++ ) {
						if (part[j] > 0.0) {
							kls[i] += contextProb * marginalDrawn[j] * Math.log(marginalDrawn[j] / marginal[j]);
						}
					}
				}

			}
		}
		

		private double getContextProbability() {
			if (this.children == null) {
				double[] vals = new double[pars.length];
				for (int i = 0; i < pars.length; i++) {
					vals[i] = pars[i].getValue() + pars[i].getLogT();
				}
				return Math.exp(Normalisation.getLogSum(vals));
			} else {
				double val = 0;
				for (int i = 0; i < children.length; i++) {
					val += children[i].getContextProbability();
				}
				return val;
			}
		}

		// TODO check next 3 methods
		private void findAndFill(double[][] fillEmptyWith, int contextLength) {
			fill(fillEmptyWith, 0, 1, contextLength);
		}

		private void fill(double[][] fillEmptyWith, int context, int power,
				int contextLength) {
			if (contextLength > 0) {
				// fill/compute the context
				contextLength--;
				for (int i = 0; i < children.length; i++) {
					children[i].fill(fillEmptyWith, context + i * power, power
							* fillEmptyWith[0].length, contextLength);
				}
			} else {
				// fill parameters
				fill(fillEmptyWith[context]);
			}
		}

		// this is the main method for filling the parameters
		private void fill(double[] distr) {
			if (children != null) {
				for (int i = 0; i < children.length; i++) {
					children[i].fill(distr);
				}
			} else {
				if (pars[pars.length - 1].isFree()) {
					if (distr.length != pars.length) {
						throw new IndexOutOfBoundsException(
								"Different number of values (" + distr.length
										+ ") than free parameters ("
										+ pars.length + ").");
					} else {
						for (int i = 0; i < pars.length; i++) {
							pars[i].setValue(Math.log(distr[i]));
						}
					}
				} else {
					for (int i = 0; i < pars.length - 1; i++) {
						pars[i].setValue(Math.log(distr[i])
								- Math.log(distr[pars.length - 1]));
					}
				}
			}
		}

		private void copy(TreeElement node) {
			if (this.children != null) {
				if (node.children != null) {
					if (this.children.length != node.children.length) {
						throw new IndexOutOfBoundsException(
								"Different number of children.");
					} else {
						for (int i = 0; i < children.length; i++) {
							this.children[i].copy(node.children[i]);
						}
					}
				} else {
					for (int i = 0; i < children.length; i++) {
						this.children[i].copy(node);
					}
				}
			} else {
				if (node.pars != null) {
					if (this.pars.length != node.pars.length) {
						throw new IndexOutOfBoundsException(
								"Different number of parameters.");
					} else {
						for (int i = 0; i < pars.length; i++) {
							pars[i].setValue(node.pars[i].getValue());
						}
					}
				} else {
					double[] vals = new double[pars.length];
					for (int i = 0; i < pars.length; i++) {
						double res = node.getLogSum(i);
						vals[i] = res;
						pars[i].setValue(res);
					}
					double norm = Normalisation.getLogSum(vals);
					for (int i = 0; i < pars.length; i++) {
						pars[i].setValue(pars[i].getValue() - norm);
					}
				}
			}
		}

		private double getLogSum(int idx) {
			if (this.children != null) {
				double[] vals = new double[this.children.length];
				for (int i = 0; i < vals.length; i++) {
					vals[i] = children[i].getLogSum(idx);
				}
				double ret = Normalisation.getLogSum(vals);
				return ret;
			} else {
				return getLogSumForLeaf(idx);
			}
		}

		private double getLogSumForLeaf(int idx) {
			return pars[idx].getValue() + pars[idx].getLogT();
		}

		private void initializeRandomly(double ess) {
			if (pars != null) {

				if (ess <= 0) {
					ess = alphabet.getAlphabetLengthAt(pars[0].getPosition());
				}
				double[] hyp = new double[pars.length];
				for (int i = 0; i < hyp.length; i++) {
					hyp[i] = ess
							/ alphabet.getAlphabetLengthAt(pars[i]
									.getPosition());
				}
				double[] temp = DirichletMRG.DEFAULT_INSTANCE.generate(
						pars.length, new DirichletMRGParams(hyp));
				for (int i = 0; i < pars.length; i++) {
					pars[i].count = temp[i];
				}
				this.normalizePlugInParameters();
				if (!pars[pars.length - 1].isFree()) {
					this.divideByUnfree();
				}
			} else {
				for (int i = 0; i < children.length; i++) {
					children[i].initializeRandomly(ess
							/ alphabet.getAlphabetLengthAt(this.contextPos));
				}
			}

		}

		private double computeGammaNorm() {
			if (children != null) {
				double val = 0;
				for (int i = 0; i < children.length; i++) {
					val += children[i].computeGammaNorm();
				}
				return val;
			} else {
				double val = 0;
				double hypSum = 0, alpha;
				for (int i = 0; i < pars.length; i++) {
					alpha = pars[i].getPseudoCount();
					hypSum += alpha;
					val -= Gamma.logOfGamma(alpha);
				}
				val += Gamma.logOfGamma(hypSum);
				return val;
			}
		}

		private double getProbFor( Sequence sequence, int offset ) {
			if(this.children != null){
				if(offset < sequence.getLength()-1){
					return children[sequence.discreteVal( offset )].getProbFor( sequence, offset+1 );
				}else{
					double val = 0.0;
					for(int i=0;i<children.length;i++){
						val += children[i].getProbFor( sequence, offset+1 );
					}
					return val;
				}
			}else{
				return getContextProbability()*pars[sequence.discreteVal( sequence.getLength()-1 )].getExpValue();
			}
		}

		private int getNumberOfParameters() {
			if(this.children != null){
				int sum = 0;
				for(int i=0;i<children.length;i++){
					sum += children[i].getNumberOfParameters();
				}
				return sum;
			}else{
				int sum = 0;
				for(int i=0;i<pars.length;i++){
					if(pars[i].isFree()){
						sum++;
					}
				}
				return sum;
			}
		}

		private int getNumberOfSamplingSteps() {
			if(this.children != null){
				int sum = 0;
				for(int i=0;i<children.length;i++){
					sum += children[i].getNumberOfSamplingSteps();
				}
				return sum;
			}else{
				return 1;
			}
		}

		private int[] getParameterIndexesForSamplingStep( int step, int offset ) {
			if(this.children != null){
				for(int i=0;i<children.length;i++){
					int currSteps = children[i].getNumberOfSamplingSteps();
					if( step < currSteps){
						return children[i].getParameterIndexesForSamplingStep( step, offset );
					}else{
						step -= currSteps;
						offset += children[i].getNumberOfParameters();
					}
				}
				return null;
			}else{
				int[] pars = new int[this.pars.length-1];
				for(int i=0;i<pars.length;i++){
					pars[i] = offset+i;
				}
				return pars;
			}
		}

	}

	/**
	 * Returns the number of sampling steps in a grouped sampling
	 * to sample all parameters of this {@link ParameterTree}.
	 * @return the number of steps
	 */
	int getNumberOfSamplingSteps(){
		return root.getNumberOfSamplingSteps();
	}

	/**
	 * Returns the number of parameters represented by this
	 * {@link ParameterTree}
	 * @return the number of parameters
	 */
	int getNumberOfParameters() {
		return root.getNumberOfParameters();
	}

	/**
	 * Returns the indexes of the parameters, incremented by <code>offset</code>, that
	 * shall be sampled in step <code>step</code> of a grouped sampling process.
	 * @param step the step
	 * @param offset the offset on the parameter indexes
	 * @return the indexes of the group of parameters
	 */
	public int[] getParameterIndexesForSamplingStep( int step, int offset ) {
		return root.getParameterIndexesForSamplingStep(step, offset);
	}

}

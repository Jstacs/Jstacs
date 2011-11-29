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
import java.util.Iterator;
import java.util.LinkedList;

import de.jstacs.InstantiableFromParameterSet;
import de.jstacs.NonParsableException;
import de.jstacs.algorithms.graphs.TopSort;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sample;
import de.jstacs.data.Sequence;
import de.jstacs.io.XMLParser;
import de.jstacs.io.ParameterSetParser.NotInstantiableException;
import de.jstacs.parameters.InstanceParameterSet;
import de.jstacs.scoringFunctions.AbstractNormalizableScoringFunction;
import de.jstacs.scoringFunctions.directedGraphicalModels.structureLearning.measures.InhomogeneousMarkov;
import de.jstacs.scoringFunctions.directedGraphicalModels.structureLearning.measures.Measure;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.IntList;

/**
 * This class implements a scoring function that is a moral directed graphical
 * model, i.e. a moral Bayesian network. This implementation also comprises well
 * known specializations of Bayesian networks like Markov models of arbitrary
 * order (including weight array matrix models (WAM) and position weight
 * matrices (PWM)) or Bayesian trees. Different structures can be achieved by
 * using the corresponding {@link Measure}, e.g. {@link InhomogeneousMarkov} for
 * Markov models of arbitrary order. <br />
 * <br />
 * 
 * This scoring function can be used in any
 * {@link de.jstacs.classifier.scoringFunctionBased.ScoreClassifier}, e.g. in a
 * {@link de.jstacs.classifier.scoringFunctionBased.msp.MSPClassifier} to learn
 * the parameters of the {@link de.jstacs.scoringFunctions.ScoringFunction}
 * using maximum conditional likelihood.
 * 
 * @author Jan Grau
 */
public class BayesianNetworkScoringFunction extends
		AbstractNormalizableScoringFunction implements
		InstantiableFromParameterSet {

	/**
	 * The parameters of the scoring function. This comprises free as well as
	 * dependent parameters.
	 */
	protected Parameter[] parameters;
	/**
	 * The trees that represent the context of the random variable (i.e.
	 * configuration of parent random variables) of the parameters.
	 */
	protected ParameterTree[] trees;
	/**
	 * Indicates if the instance has been trained.
	 */
	protected boolean isTrained;
	/**
	 * The equivalent sample size.
	 */
	protected double ess;
	/**
	 * The number of free parameters.
	 */
	protected Integer numFreePars;
	/**
	 * Used internally. Mapping from indexes of free parameters to indexes of
	 * all parameters.
	 */
	protected int[] nums;
	/**
	 * {@link Measure} that defines the network structure.
	 */
	protected Measure structureMeasure;
	/**
	 * Indicates if plug-in parameters, i.e. generative (MAP) parameters shall
	 * be used upon initialization.
	 */
	protected boolean plugInParameters;
	/**
	 * The network structure, used internally.
	 */
	protected int[][] order;
	/**
	 * Normalization constant to obtain normalized probabilities.
	 */
	protected Double logNormalizationConstant;
	/**
	 * The roots of the network structure, i.e. random variables without
	 * parents.
	 */
	private int[] roots;
	/**
	 * Indicates if all or only the free parameters should be used for
	 * optimization.
	 */
	private boolean freeParams;
	/**
	 * Precomputed Gamma-normalization.
	 */
	private Double gammaNorm;

	private BayesianNetworkScoringFunctionParameterSet parameterSet;

	/**
	 * Creates a new {@link BayesianNetworkScoringFunction} that has neither
	 * been initialized nor trained.
	 * 
	 * @param alphabet
	 *            the alphabet of the scoring function boxed in an
	 *            {@link AlphabetContainer}, e.g
	 *            <code>new AlphabetContainer(new DNAAlphabet())</code>
	 * @param length
	 *            the length of the scoring function, i.e. the length of the
	 *            sequences this scoring function can handle
	 * @param ess
	 *            the equivalent sample size
	 * @param plugInParameters
	 *            indicates if plug-in parameters, i.e. generative (MAP)
	 *            parameters, shall be used upon initialization
	 * @param structureMeasure
	 *            the {@link Measure} used for the structure, e.g.
	 *            {@link InhomogeneousMarkov}
	 * @throws Exception
	 *             if the length of the scoring function is not admissible (<=0)
	 *             or the alphabet is not discrete
	 */
	public BayesianNetworkScoringFunction(AlphabetContainer alphabet,
			int length, double ess, boolean plugInParameters,
			Measure structureMeasure) throws Exception {
		super(alphabet, length);
		if (!alphabet.isDiscrete()) {
			throw new Exception("Only defined on discrete alphabets.");
		}
		if (length <= 0) {
			throw new Exception("Inconsistent length (" + length + ").");
		}
		isTrained = false;
		this.ess = ess;
		this.plugInParameters = plugInParameters;
		this.structureMeasure = structureMeasure;
		this.logNormalizationConstant = null;
	}

	/**
	 * Creates a new {@link BayesianNetworkScoringFunction} that has neither
	 * been initialized nor trained from a
	 * {@link BayesianNetworkScoringFunctionParameterSet}.
	 * 
	 * @param parameters
	 *            the parameter set
	 * 
	 * @throws NotInstantiableException
	 *             if the {@link BayesianNetworkScoringFunction} could not be
	 *             instantiated from the
	 *             {@link BayesianNetworkScoringFunctionParameterSet}
	 * @throws Exception
	 *             if the length of the scoring function is not admissible (<=0)
	 *             or the alphabet is not discrete
	 */
	public BayesianNetworkScoringFunction(
			BayesianNetworkScoringFunctionParameterSet parameters)
			throws NotInstantiableException, Exception {
		this(parameters.getAlphabetContainer(), parameters.getLength(), parameters
				.getEss(), parameters.getPlugInParameters(), parameters
				.getMeasure());
		this.parameterSet = parameters;
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Recreates a {@link BayesianNetworkScoringFunction} from its XML
	 * representation as saved by the method {@link #toXML()}.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML code could not be parsed
	 */
	public BayesianNetworkScoringFunction(StringBuffer xml)
			throws NonParsableException {
		super(xml);
		this.logNormalizationConstant = null;
		this.gammaNorm = null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.scoringFunctions.AbstractNormalizableScoringFunction#clone()
	 */
	@SuppressWarnings("unchecked")
	@Override
	public BayesianNetworkScoringFunction clone()
			throws CloneNotSupportedException {
		BayesianNetworkScoringFunction clone = (BayesianNetworkScoringFunction) super
				.clone();
		if (trees != null) {
			clone.trees = new ParameterTree[trees.length];
			for (int i = 0; i < trees.length; i++) {
				clone.trees[i] = trees[i].clone();
			}

			LinkedList<Parameter>[] parTemp = new LinkedList[trees.length];
			int num = 0;
			for (int i = 0; i < clone.trees.length; i++) {
				parTemp[i] = clone.trees[i].linearizeParameters();
				num += parTemp[i].size();
			}
			clone.parameters = new Parameter[num];
			num = 0;
			Iterator<Parameter> it = null;
			for (int i = 0; i < parTemp.length; i++) {
				it = parTemp[i].iterator();
				while (it.hasNext()) {
					clone.parameters[num++] = it.next();
				}
			}
			clone.nums = nums.clone();
		} else {
			clone.trees = null;
			clone.nums = null;
			clone.parameters = null;
		}
		clone.structureMeasure = structureMeasure.clone();
		clone.logNormalizationConstant = null;
		return clone;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @seede.jstacs.scoringFunctions.NormalizableScoringFunction#
	 * getPartialNormalizationConstant(int)
	 */
	public double getLogPartialNormalizationConstant(int parameterIndex)
			throws Exception {
		if (logNormalizationConstant == null) {
			precomputeNormalization();
		}
		if (parameterIndex < nums.length) {
			Parameter p = parameters[nums[parameterIndex]];
			int pos = p.getPosition();

			boolean[] notRoot = new boolean[trees.length];

			for (int i = 0; i < trees.length; i++) {
				if (trees[i].getNumberOfParents() > 0) {
					notRoot[i] = true;
				}
			}

			double val = p.getLogPartialNormalizer();
			for (int i = 0; i < notRoot.length; i++) {
				if (!notRoot[i] && i != roots[pos]) {
					val += trees[i].forward(trees);
				}
			}
			return val;

		} else {
			throw new Exception("Parameter index out of bounds");
		}
	}
	

	@SuppressWarnings("unchecked")
	private int[][] getFirstChildrenAndFirstParents(int[][] parents)
			throws Exception {
		LinkedList<Integer>[] fc = new LinkedList[parents.length];
		boolean test = false;

		for (int i = 0; i < fc.length; fc[i++] = new LinkedList<Integer>())
			;

		int[][] erg = new int[parents.length + 1][];
		erg[parents.length] = new int[parents.length];
		Arrays.fill(erg[parents.length], -1);

		for (int i = 0; i < parents.length; i++) {
			// no parents, always OK
			test = parents[i].length < 2;
			for (int j = 0; j < parents[i].length - 1; j++) {
				if (testInclude(parents[i], parents[parents[i][j]])) {
					fc[parents[i][j]].add(i);
					erg[parents.length][i] = parents[i][j];
					test = true;
				}
			}
			if (!test) {
				throw new Exception("Structure is no moral graph!");
			}
		}

		for (int i = 0; i < fc.length; i++) {
			erg[i] = new int[fc[i].size()];
			for (int j = 0; j < erg[i].length; j++) {
				erg[i][j] = fc[i].poll();
			}
		}
		return erg;
	}

	private boolean testInclude(int[] parentsOfChild, int[] parentsOfParent) {

		for (int i = 0; i < parentsOfChild.length - 1; i++) {
			boolean found = false;
			innerloop: for (int j = 0; j < parentsOfParent.length; j++) {
				if (parentsOfChild[i] == parentsOfParent[j]) {
					found = true;
					break innerloop;
				}
			}
			if (!found) {
				return false;
			}
		}
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.scoringFunctions.ScoringFunction#initializeFunction(int,
	 * boolean, de.jstacs.data.Sample[], double[][])
	 */
	public void initializeFunction(int index, boolean freeParams,
			Sample[] data, double[][] weights) throws Exception {
		if (data[index] != null && data[index].getElementLength() != length) {
			throw new Exception("Data has wrong length.");
		}
		this.freeParams = freeParams;
		if (weights == null) {
			weights = new double[data.length][];
			for (int i = 0; i < data.length; i++) {
				weights[i] = new double[data[i].getNumberOfElements()];
				Arrays.fill(weights[i], 1.0);
			}
		}
		Sample[] data2 = data;
		double[][] weights2 = weights;
		if (data.length != 2) {
			data2 = new Sample[2];
			weights2 = new double[2][];
			data2[0] = data[index];
			weights2[0] = weights[index];
			boolean[] in = new boolean[data.length];
			Arrays.fill(in, true);
			in[index] = false;

			data2[1] = Sample.union(data, in);
			if (data2[1] != null) {
				weights2[1] = new double[data2[1].getNumberOfElements()];
				int off = 0;
				for (int i = 0; i < weights.length; i++) {
					if (in[i]) {
						System.arraycopy(weights[i], 0, weights2[1], off,
								weights[i].length);
						off += weights[i].length;
					}
				}
			} else {
				weights2[1] = null;
			}
		}

		createTrees(data2, weights2);

		if (plugInParameters) {
			setPlugInParameters(index, freeParams, data, weights);
		} else {
			for (int i = 0; i < parameters.length; i++) {
				if(freeParams){
					parameters[i].setValue( 0 );
				}else{
					parameters[i].setValue(-Math.log(alphabets
						.getAlphabetLengthAt(parameters[i].position)));
				}
			}
		}

		isTrained = true;
		logNormalizationConstant = null;
		// System.out.println(this.toString());
	}

	/**
	 * Creates the tree structures that represent the context (array
	 * {@link #trees}) and the parameter objects {@link #parameters} using the
	 * given {@link Measure} {@link #structureMeasure}.
	 * 
	 * @param data2
	 *            the data that is used to compute the structure
	 * @param weights2
	 *            the weights on the sequences in <code>data2</code>
	 * 
	 * @throws Exception
	 *             if the structure is no moral graph or if the lengths of data
	 *             and scoring function do not match or other problems
	 *             concerning the data occur
	 */
	protected void createTrees(Sample[] data2, double[][] weights2)
			throws Exception {
		int[][] parents = structureMeasure.getParents(data2[0], data2[1],
				weights2[0], weights2[1], this.getLength());

		this.order = TopSort.getTopologicalOrder(parents);

		int[][] firstChildrenAndFirstParents = getFirstChildrenAndFirstParents(parents);

		numFreePars = 0;
		int numPars = 0;
		int[] numContextsPos = new int[this.getLength()];
		trees = new ParameterTree[this.getLength()];
		int numContexts = 0;

		int[] contextPoss = null;

		for (int i = 0; i < parents.length; i++) {
			numContextsPos[i] = 1;
			for (int j = 0; j < parents[i].length - 1; j++) {
				numContextsPos[i] *= alphabets
						.getAlphabetLengthAt(parents[i][j]);
			}
			numFreePars += numContextsPos[i]
					* ((int) getAlphabetContainer().getAlphabetLengthAt(i) - 1);
			numPars += numContextsPos[i]
					* getAlphabetContainer().getAlphabetLengthAt(i);
			numContexts += numContextsPos[i];
			contextPoss = new int[parents[i].length - 1];
			for (int j = 0; j < contextPoss.length; j++) {
				contextPoss[j] = parents[i][parents[i].length - j - 2];
			}
			trees[i] = new ParameterTree(i, contextPoss,
					getAlphabetContainer(),
					firstChildrenAndFirstParents[parents.length][i],
					firstChildrenAndFirstParents[i]);
		}

		this.parameters = new Parameter[numPars];
		if (!freeParams) {
			numFreePars = numPars;
		}
		this.nums = new int[numFreePars];

		int curr = 0;
		int free = 0;
		// [position][contextNum][depth][0: pos, 1-...: symbols]
		int[][][][] contexts = new int[this.getLength()][][][];

		for (int i = 0; i < parents.length; i++) {
			contexts[i] = new int[numContextsPos[i]][parents[i].length - 1][];

			fillContexts(0, contexts[i], 0, parents[i]);

			for (int j = 0; j < contexts[i].length; j++) {
				int all = 1;
				int act = 1;

				for (int k = 0; k < contexts[i][j].length; k++) {
					all *= getAlphabetContainer().getAlphabetLengthAt(
							contexts[i][j][k][0]);
					act *= contexts[i][j][k].length - 1;
				}
				for (byte a = 0; a < getAlphabetContainer()
						.getAlphabetLengthAt(i); a++) {
					if (a < getAlphabetContainer().getAlphabetLengthAt(i) - 1
							|| (!freeParams)) {
						parameters[curr] = new Parameter(
								free,
								a,
								i,
								contexts[i][j],
								(double) (act * this.ess)
										/ (double) (all * getAlphabetContainer()
												.getAlphabetLengthAt(i)), true);
					} else {
						parameters[curr] = new Parameter(
								-1,
								a,
								i,
								contexts[i][j],
								(double) (act * this.ess)
										/ (double) (all * getAlphabetContainer()
												.getAlphabetLengthAt(i)), false);
					}
					trees[i].setParameterFor(a, contexts[i][j],
							parameters[curr]);
					if (parameters[curr].isFree()) {
						nums[free++] = curr;
					}
					curr++;
				}
			}
		}

		roots = new int[trees.length];
		for (int i = 0; i < trees.length; i++) {
			int fp = i;
			while (trees[fp].getFirstParent() != -1) {
				fp = trees[fp].getFirstParent();
			}
			roots[i] = fp;
		}
		logNormalizationConstant = null;
		gammaNorm = null;
	}

	/**
	 * Computes and sets the plug-in parameters (MAP estimated parameters) from
	 * <code>data</code> using <code>weights</code>.
	 * 
	 * @param index
	 *            the index of the class the scoring function is responsible
	 *            for, the parameters are estimated from
	 *            <code>data[index]</code> and <code>weights[index]</code>
	 * @param freeParameters
	 *            indicates if only the free parameters or all parameters should
	 *            be used, this also affects the initialization
	 * @param data
	 *            the data used for initialization
	 * @param weights
	 *            the weights on the data
	 */
	protected void setPlugInParameters(int index, boolean freeParameters,
			Sample[] data, double[][] weights) {
		if (data[index] != null) {
			for (int i = 0; i < data[index].getNumberOfElements(); i++) {
				for (int j = 0; j < trees.length; j++) {
					trees[j].addCount(data[index].getElementAt(i), 0,
							weights[index][i]);
				}
			}

			for (int i = 0; i < trees.length; i++) {
				trees[i].normalizePlugInParameters();
				if (freeParameters) {
					trees[i].divideByUnfree();
				}
			}
		}
	}

	private int fillContexts(int offset, int[][][] contexts, int depth,
			int[] parents) {
		int tempOffset = offset;

		if (depth < parents.length - 1) {
			for (int i = 0; i < alphabets
					.getAlphabetLengthAt(parents[parents.length - depth - 2]); i++) {
				offset = tempOffset;
				tempOffset = fillContexts(offset, contexts, depth + 1, parents);
				for (int j = offset; j < tempOffset; j++) {
					contexts[j][depth] = new int[] {
							parents[parents.length - depth - 2], i };
				}
			}
			return tempOffset;
		} else {
			return offset + 1;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.scoringFunctions.AbstractNormalizableScoringFunction#fromXML
	 * (java.lang.StringBuffer)
	 */
	@SuppressWarnings("unchecked")
	@Override
	protected void fromXML(StringBuffer source) throws NonParsableException {
		source = XMLParser.extractForTag(source, "bayesianNetworkSF");
		alphabets = XMLParser.extractObjectForTags(source, "alphabets", AlphabetContainer.class);
		length = XMLParser.extractObjectForTags(source, "length", int.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
		trees = XMLParser.extractObjectForTags(source, "trees", ParameterTree[].class );
		if (trees.length == 0) {
			trees = null;
			this.parameters = null;
		} else {
			LinkedList<Parameter>[] parTemp = new LinkedList[trees.length];
			int num = 0;
			for (int i = 0; i < trees.length; i++) {
				trees[i].setAlphabet( alphabets );
				parTemp[i] = trees[i].linearizeParameters();
				num += parTemp[i].size();
			}
			this.parameters = new Parameter[num];
			num = 0;
			Iterator<Parameter> it = null;
			for (int i = 0; i < parTemp.length; i++) {
				it = parTemp[i].iterator();
				while (it.hasNext()) {
					parameters[num++] = it.next();
				}
			}
		}
		isTrained = XMLParser.extractObjectForTags(source, "isTrained", boolean.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
		ess = XMLParser.extractObjectForTags(source, "ess", double.class );

		numFreePars = XMLParser.extractObjectForTags(source, "numFreePars", Integer.class );
		/*
		String numFreeParsTemp = XMLParser.extractObjectForTags(source, "numFreePars", String.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
		if (numFreeParsTemp.equals("null")) {
			numFreePars = null;
		} else {
			numFreePars = Integer.parseInt(numFreeParsTemp);
		}
		*/
		nums = XMLParser.extractObjectForTags(source, "nums", int[].class );
		if (nums.length == 0) {
			nums = null;
		}
		structureMeasure = XMLParser.extractObjectForTags(source, "structureMeasure", Measure.class );
		order = XMLParser.extractObjectForTags(source, "order", int[][].class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
		if (order.length == 0) {
			order = null;
		}
		plugInParameters = XMLParser.extractObjectForTags(source, "plugInParameters", boolean.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
		roots = XMLParser.extractObjectForTags(source, "roots", int[].class );
		if (roots.length == 0) {
			roots = null;
		}
		freeParams = XMLParser.extractObjectForTags(source, "freeParams", boolean.class );
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		if (trees != null) {
			if (logNormalizationConstant == null) {
				precomputeNormalization();
			}
			StringBuffer buf = new StringBuffer();
			for (int i = 0; i < trees.length; i++) {
				buf.append(trees[i].toString());
				buf.append("\n");
			}
			return buf.toString();
		} else {
			return "BNScoringFunction of length " + length
					+ ": not initialized";
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.scoringFunctions.ScoringFunction#getInstanceName()
	 */
	public String getInstanceName() {
		return structureMeasure.getInstanceName();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.scoringFunctions.ScoringFunction#getLogScore(de.jstacs.data
	 * .Sequence, int)
	 */
	public double getLogScoreFor(Sequence seq, int start) {
		double prob = 0;
		for (int i = 0; i < trees.length; i++) {
			prob += trees[i].getParameterFor(seq, start).getValue();
		}
		return prob;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.scoringFunctions.ScoringFunction#getLogScoreAndPartialDerivation
	 * (de.jstacs.data.Sequence, int, de.jstacs.utils.IntList,
	 * de.jstacs.utils.DoubleList)
	 */
	public double getLogScoreAndPartialDerivation(Sequence seq, int start,
			IntList indices, DoubleList partialDer) {
		double logScore = 0;
		for (int i = 0; i < trees.length; i++) {
			Parameter par = trees[i].getParameterFor(seq, start);
			if (par.isFree()) {
				indices.add(par.getIndex());
				partialDer.add(1);
			}
			logScore += par.getValue();
		}
		return logScore;
	}

	
	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.scoringFunctions.NormalizableScoringFunction#
	 * getNormalizationConstant()
	 */
	public double getLogNormalizationConstant() throws RuntimeException {
		if (logNormalizationConstant == null) {
			precomputeNormalization();
		}
		return logNormalizationConstant;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.scoringFunctions.ScoringFunction#getNumberOfParameters()
	 */
	public int getNumberOfParameters() {
		if (nums == null) {
			return UNKNOWN;
		} else {
			return nums.length;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.scoringFunctions.ScoringFunction#setParameters(double[],
	 * int)
	 */
	public void setParameters(double[] params, int start) {
		for (int i = 0; i < nums.length; i++) {
			this.parameters[nums[i]].setValue(params[i + start]);
		}
		this.logNormalizationConstant = null;
	}

	/**
	 * Pre-computes all normalization constants.
	 */
	protected void precomputeNormalization() {

		boolean[] notRoot = new boolean[trees.length];

		for (int i = 0; i < trees.length; i++) {
			trees[i].invalidateNormalizers();
		}
		for (int i = 0; i < trees.length; i++) {
			if (trees[i].getNumberOfParents() > 0) {
				notRoot[i] = true;
			}
			if (trees[i].isLeaf()) {
				trees[i].backward(trees, order);
			}
		}
		double val = 0;
		for (int i = 0; i < notRoot.length; i++) {
			if (!notRoot[i]) {
				val += trees[i].forward(trees);
			}
		}
		logNormalizationConstant = val;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.scoringFunctions.ScoringFunction#getCurrentParameterValues()
	 */
	public double[] getCurrentParameterValues() throws Exception {
		double[] pars = new double[nums.length];
		for (int i = 0; i < pars.length; i++) {
			pars[i] = parameters[nums[i]].getValue();
		}
		return pars;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.Storable#toXML()
	 */
	public StringBuffer toXML() {
		StringBuffer source = new StringBuffer();
		XMLParser.appendObjectWithTags(source, alphabets, "alphabets");
		XMLParser.appendObjectWithTags(source, length, "length");
		XMLParser.appendObjectWithTags(source, trees, "trees");
		
		XMLParser.appendObjectWithTags(source, isTrained, "isTrained");
		XMLParser.appendObjectWithTags(source, ess, "ess");
		XMLParser.appendObjectWithTags(source, numFreePars, "numFreePars");
		/*
		if (numFreePars == null) {
			XMLParser.appendObjectWithTags(source, "null", "numFreePars");
		} else {
			XMLParser.appendObjectWithTags(source, numFreePars, "numFreePars");
		}
		*/
		if (nums == null) {
			XMLParser.appendObjectWithTags(source, new int[0], "nums");
		} else {
			XMLParser.appendObjectWithTags(source, nums, "nums");
		}
		XMLParser.appendObjectWithTags(source, structureMeasure, "structureMeasure");
		XMLParser.appendObjectWithTags(source, plugInParameters, "plugInParameters");
		if (order == null) {
			XMLParser.appendObjectWithTags(source, new int[0][0], "order");
		} else {
			XMLParser.appendObjectWithTags(source, order, "order");
		}
		if (roots == null) {
			XMLParser.appendObjectWithTags(source, new int[0], "roots");
		} else {
			XMLParser.appendObjectWithTags(source, roots, "roots");
		}
		XMLParser.appendObjectWithTags(source, freeParams, "freeParams");
		XMLParser.addTags(source, "bayesianNetworkSF");
		return source;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.scoringFunctions.NormalizableScoringFunction#getLogPriorTerm()
	 */
	public double getLogPriorTerm() {
		if (gammaNorm == null) {
			computeGammaNorm();
		}

		double val = 0;
		for (int i = 0; i < nums.length; i++) {
			val += this.parameters[nums[i]].getValue()
					* parameters[nums[i]].getPseudoCount();
		}
		return val + gammaNorm;
	}

	private void computeGammaNorm() {
		this.gammaNorm = 0.0;
		for (int i = 0; i < trees.length; i++) {
			gammaNorm += trees[i].computeGammaNorm();
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @seede.jstacs.scoringFunctions.NormalizableScoringFunction#
	 * addGradientOfLogPriorTerm(double[], int)
	 */
	public void addGradientOfLogPriorTerm(double[] grad, int start) {
		for (int i = 0; i < nums.length; i++) {
			grad[i + start] += this.parameters[nums[i]].getPseudoCount();
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.scoringFunctions.NormalizableScoringFunction#getEss()
	 */
	public double getESS() {
		return ess;
	}

	/**
	 * Returns the position in the sequence the parameter <code>index</code> is
	 * responsible for.
	 * 
	 * @param index
	 *            the index of the parameter
	 * 
	 * @return the position in the sequence
	 */
	public int getPositionForParameter(int index) {
		return parameters[nums[index]].getPosition();
	}
	
	/**
	 * Returns the probability of <code>kmer</code> for all possible positions in this {@link BayesianNetworkScoringFunction} starting at position <code>kmer.getLength()-1<code>.
	 * @param kmer the k-mer
	 * @return the position-dependent probabilities of this k-mer for position <code>kmer.getLength()-1<code> to <code>{@link BayesianNetworkScoringFunction#getLength()}-1</code>
	 * @throws Exception if the method is called for non-Markov model structures
	 */
	public double[] getPositionDependentKMerProb(Sequence kmer) throws Exception{
		if(!(structureMeasure instanceof InhomogeneousMarkov)){
			throw new Exception("Only implemented for IMMs");
		}
		precomputeNormalization();
		for(int i=0;i<trees.length;i++){
			trees[i].normalizeParameters();
		}
		precomputeNormalization();
		
		double[] prof = new double[trees.length-kmer.getLength()+1];
		
		
		Arrays.fill( prof, 1 );
		
		for(int i=0;i<trees.length-kmer.getLength()+1;i++){
			prof[i] *= trees[i+kmer.getLength()-1].getProbFor(kmer);
		}
		return prof;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @seede.jstacs.scoringFunctions.NormalizableScoringFunction#
	 * getSizeOfEventSpaceForRandomVariablesOfParameter(int)
	 */
	public int getSizeOfEventSpaceForRandomVariablesOfParameter(int index) {
		Parameter par = parameters[nums[index]];
		int size = (int) alphabets.getAlphabetLengthAt(par.getPosition());
		int[][] cont = par.context;
		for (int i = 0; i < cont.length; i++) {
			size *= alphabets.getAlphabetLengthAt(cont[i][0]);
		}
		return size;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.scoringFunctions.ScoringFunction#initializeFunctionRandomly
	 * (boolean)
	 */
	public void initializeFunctionRandomly(boolean freeParams) throws Exception {
		if (!(structureMeasure instanceof InhomogeneousMarkov)) {
			throw new Exception("Not implemented");
		}
		this.freeParams = freeParams;
		boolean temp = this.plugInParameters;
		this.plugInParameters = false;
		this.initializeFunction(0, freeParams, new Sample[] { null, null },
				new double[][] { null, null });
		this.plugInParameters = temp;

		for (int i = 0; i < trees.length; i++) {
			trees[i].initializeRandomly(ess);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.scoringFunctions.ScoringFunction#isInitialized()
	 */
	public boolean isInitialized() {
		return trees != null;
	}

	/**
	 * If this {@link BayesianNetworkScoringFunction} is a PWM, i.e.
	 * {@link #structureMeasure}=new {@link InhomogeneousMarkov}(0)}}, this
	 * method returns the normalized PWM as a <code>double</code> array of
	 * dimension {@link #getLength()} x size-of-alphabet.
	 * 
	 * @return the PWM as a two-dimensional array
	 * 
	 * @throws Exception
	 *             if this method is called for a
	 *             {@link BayesianNetworkScoringFunction} that is not a PWM
	 */
	public double[][] getPWM() throws Exception {
		if (!(this.structureMeasure instanceof InhomogeneousMarkov)) {
			throw new Exception();
		}
		if (logNormalizationConstant == null) {
			precomputeNormalization();
		}
		double[][] pwm = new double[trees.length][];
		for (int i = 0; i < trees.length; i++) {
			pwm[i] = new double[(int) alphabets.getAlphabetLengthAt(i)];
			trees[i].insertProbs(pwm[i]);
		}
		return pwm;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.InstantiableFromParameterSet#getCurrentParameterSet()
	 */
	public InstanceParameterSet getCurrentParameterSet() throws Exception {
		if (this.parameterSet != null) {
			return this.parameterSet;
		} else {
			return new BayesianNetworkScoringFunctionParameterSet(
					getAlphabetContainer(), getLength(), ess, plugInParameters,
					structureMeasure);
		}
	}
}

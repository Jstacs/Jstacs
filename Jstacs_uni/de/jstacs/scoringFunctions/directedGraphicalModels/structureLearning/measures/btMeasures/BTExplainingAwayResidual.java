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

package de.jstacs.scoringFunctions.directedGraphicalModels.structureLearning.measures.btMeasures;

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.algorithms.graphs.MST;
import de.jstacs.data.Sample;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.InstanceParameterSet;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.scoringFunctions.directedGraphicalModels.structureLearning.measures.Measure;

/**
 * Structure learning {@link Measure} that computes a maximum spanning tree
 * based on the explaining away residual and uses the resulting tree structure
 * as structure of a Bayesian tree (special case of a Bayesian network) in a
 * {@link de.jstacs.scoringFunctions.directedGraphicalModels.BayesianNetworkScoringFunction}
 * .
 * 
 * @author Jan Grau
 */
public class BTExplainingAwayResidual extends Measure {

	private BTExplainingAwayResidualParameterSet parameters;
	private double[] ess;

	/**
	 * Creates a new explaining away residual Bayesian tree {@link Measure}.
	 * 
	 * @param ess
	 *            the equivalent sample sizes (ess) for the classes
	 */
	public BTExplainingAwayResidual(double[] ess) {
		this.ess = ess;
	}

	/**
	 * Creates a new {@link BTExplainingAwayResidual} from the corresponding
	 * {@link InstanceParameterSet} <code>parameters</code>.
	 * 
	 * @param parameters
	 *            the corresponding parameters
	 */
	public BTExplainingAwayResidual(
			BTExplainingAwayResidualParameterSet parameters) {
		this(parameters.getEss());
		this.parameters = parameters;
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Recreates a {@link BTExplainingAwayResidual} from is XML representation
	 * as returned by {@link #toXML()}.
	 * 
	 * @param buf
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML code could not be parsed
	 */
	public BTExplainingAwayResidual(StringBuffer buf)
			throws NonParsableException {
		buf = XMLParser.extractForTag(buf, "btExplainingAwayResidual");
		ess = XMLParser.extractObjectForTags(buf, "ess", double[].class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.scoringFunctions.directedGraphicalModels.structureLearning.
	 * measures.Measure#clone()
	 */
	@Override
	public BTExplainingAwayResidual clone() throws CloneNotSupportedException {
		BTExplainingAwayResidual clone = (BTExplainingAwayResidual) super
				.clone();
		clone.ess = ess.clone();
		return clone;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.scoringFunctions.directedGraphicalModels.structureLearning.
	 * measures.Measure#getInstanceName()
	 */
	@Override
	public String getInstanceName() {
		return "Bayesian tree with explaining away residual";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.scoringFunctions.directedGraphicalModels.structureLearning.
	 * measures.Measure#getParents(de.jstacs.data.Sample, de.jstacs.data.Sample,
	 * double[], double[], int)
	 */
	@Override
	public int[][] getParents(Sample fg, Sample bg, double[] weightsFg,
			double[] weightsBg, int length) throws Exception {

		double[][][][] statFg = getStatistics(fg, weightsFg, length, ess[0]);
		double[][][][] statBg = getStatistics(bg, weightsBg, length, ess[1]);
		double nFg = sum(weightsFg) + ess[0];
		double nBg = sum(weightsBg) + ess[1];

		double[][] ear = getEAR(statFg, statBg, nFg, nBg);

		int[][] p = MST.kruskal(ear);

		int[][] parents = new int[length][1];
		for (int i = 0; i < parents.length; i++) {
			parents[i][0] = i;
		}
		for (int i = 0; i < p.length; i++) {
			int idx = p[i][1];
			parents[idx] = new int[2];
			parents[idx][0] = p[i][0];
			parents[idx][1] = idx;

		}
		return parents;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.Storable#toXML()
	 */
	public StringBuffer toXML() {
		StringBuffer buf = new StringBuffer();
		XMLParser.appendObjectWithTags(buf, ess, "ess");
		XMLParser.addTags(buf, "btExplainingAwayResidual");
		return buf;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.InstantiableFromParameterSet#getCurrentParameterSet()
	 */
	public InstanceParameterSet getCurrentParameterSet() throws Exception {
		if (parameters != null) {
			return parameters;
		} else {
			return new BTExplainingAwayResidualParameterSet(ess);
		}
	}

	/**
	 * Class for the parameters of a {@link BTExplainingAwayResidual} structure
	 * {@link Measure}.
	 * 
	 * @author Jan Grau
	 */
	public static class BTExplainingAwayResidualParameterSet extends
			InstanceParameterSet {

		/**
		 * Creates a new {@link BTExplainingAwayResidualParameterSet} with empty
		 * parameter values.
		 */
		public BTExplainingAwayResidualParameterSet() {
			super(BTExplainingAwayResidual.class);
		}

		/**
		 * Creates a new {@link BTExplainingAwayResidualParameterSet} with the
		 * parameter for the equivalent sample sizes set to <code>ess</code>.
		 * 
		 * @param ess
		 *            the equivalent sample sizes for the foreground class and
		 *            the background, i.e. the background class or (in case of
		 *            more than two classes) all non-foreground classes
		 * 
		 * @throws Exception
		 *             if the parameters could not be created or set
		 */
		public BTExplainingAwayResidualParameterSet(double[] ess)
				throws Exception {
			super(BTExplainingAwayResidual.class);
			loadParameters();
			parameters.get(0).setValue(ess[0]);
			parameters.get(1).setValue(ess[1]);
		}

		/**
		 * The standard constructor for the interface {@link de.jstacs.Storable}
		 * . Creates a new {@link BTExplainingAwayResidualParameterSet} from its
		 * XML representation as defined by the {@link de.jstacs.Storable}
		 * interface.
		 * 
		 * @param representation
		 *            the XML code as {@link StringBuffer}
		 * 
		 * @throws NonParsableException
		 *             if the XML representation could not be parsed
		 */
		public BTExplainingAwayResidualParameterSet(StringBuffer representation)
				throws NonParsableException {
			super(representation);
		}

		/**
		 * Returns the equivalent sample sizes (ess) defined by this set of
		 * parameters.
		 * 
		 * @return the equivalent sample sizes (ess)
		 */
		public double[] getEss() {
			return new double[] { (Double) parameters.get(0).getValue(),
					(Double) parameters.get(1).getValue() };
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see de.jstacs.parameters.InstanceParameterSet#getInstanceComment()
		 */
		@Override
		public String getInstanceComment() {
			return "Bayesian tree - explaining away residual";
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see de.jstacs.parameters.InstanceParameterSet#getInstanceName()
		 */
		@Override
		public String getInstanceName() {
			return "Build a Bayesian tree using explaining away residual as structure measure.";
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see de.jstacs.parameters.ParameterSet#loadParameters()
		 */
		@Override
		protected void loadParameters() throws Exception {
			initParameterList();
			this.parameters.add(new SimpleParameter(DataType.DOUBLE,
					"Foreground ESS",
					"The equivalent sample size for the foreground.", true));
			this.parameters
					.add(new SimpleParameter(
							DataType.DOUBLE,
							"Background ESS",
							"The equivalent sample size for the background,"
									+ " i.e. the background class or (in case of more than two classes) all non-foreground classes.",
							true));
		}

	}

}

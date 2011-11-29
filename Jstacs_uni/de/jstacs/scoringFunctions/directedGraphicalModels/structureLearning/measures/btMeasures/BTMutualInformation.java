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
import de.jstacs.parameters.EnumParameter;
import de.jstacs.parameters.InstanceParameterSet;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.scoringFunctions.directedGraphicalModels.structureLearning.measures.Measure;

/**
 * Structure learning {@link Measure} that computes a maximum spanning tree
 * based on mutual information and uses the resulting tree structure as
 * structure of a Bayesian tree (special case of a Bayesian network) in a
 * {@link de.jstacs.scoringFunctions.directedGraphicalModels.BayesianNetworkScoringFunction}
 * .
 * 
 * @author Jan Grau
 */
public class BTMutualInformation extends Measure {

	private BTMutualInformationParameterSet parameters;
	private DataSource clazz;
	private double[] ess;

	/**
	 * {@link Enum} defining the possible sources of data to compute the mutual
	 * information.
	 * 
	 * @author Jan Grau
	 * 
	 */
	public static enum DataSource {
		/**
		 * Compute mutual information only from foreground data.
		 */
		FG,
		/**
		 * Compute mutual information only from background data.
		 */
		BG,
		/**
		 * Use both data sets (forground and background data) to compute the
		 * mutual information.
		 */
		BOTH
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Recreates a {@link BTMutualInformation} from is XML representation as
	 * returned by {@link #toXML()}.
	 * 
	 * @param buf
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML code could not be parsed
	 */
	public BTMutualInformation(StringBuffer buf) throws NonParsableException {
		buf = XMLParser.extractForTag(buf, "btMutualInformation");
		clazz = XMLParser.extractObjectForTags(buf, "clazz", DataSource.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
		ess = XMLParser.extractObjectForTags(buf, "ess", double[].class );
	}

	/**
	 * Creates a new mutual information Bayesian tree {@link Measure}.
	 * 
	 * @param clazz
	 *            the classes used for computation of mutual information, as
	 *            defined by {@link DataSource}
	 * @param ess
	 *            the equivalent sample sizes for both classes
	 */
	public BTMutualInformation(DataSource clazz, double[] ess) {
		this.clazz = clazz;
		this.ess = ess.clone();
	}

	/**
	 * Creates a new {@link BTMutualInformation} from the corresponding
	 * {@link InstanceParameterSet} <code>parameters</code>.
	 * 
	 * @param parameters
	 *            the corresponding parameters
	 */
	public BTMutualInformation(BTMutualInformationParameterSet parameters) {
		this(parameters.getClazz(), parameters.getEss());
		this.parameters = parameters;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.scoringFunctions.directedGraphicalModels.structureLearning.
	 * measures.Measure#clone()
	 */
	@Override
	public BTMutualInformation clone() throws CloneNotSupportedException {
		BTMutualInformation clone = (BTMutualInformation) super.clone();
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
		String str = "Bayesian tree with mutual information of";
		if (clazz == DataSource.FG) {
			return str + " foreground";
		} else if (clazz == DataSource.BG) {
			return str + " background";
		} else {
			return str + " foreground and background";
		}
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
		Sample data = null;
		double[] weights = null;
		double ess2 = 0;
		if (clazz == DataSource.FG) {
			data = fg;
			weights = weightsFg;
			ess2 = ess[0];
		} else if (clazz == DataSource.BG) {
			data = bg;
			weights = weightsBg;
			ess2 = ess[1];
		} else {
			data = Sample.union(fg, bg);
			weights = union(new double[][] { weightsFg, weightsBg });
			ess2 = ess[0] + ess[1];
		}
		double[][][][] stat = getStatistics(data, weights, length, ess2);
		double[][] mi = getMI(stat, sum(weights) + ess2);

		int[][] p = null;

		p = MST.kruskal(mi);

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
		XMLParser.appendObjectWithTags(buf, clazz, "clazz");
		XMLParser.appendObjectWithTags(buf, ess, "ess");
		XMLParser.addTags(buf, "btMutualInformation");
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
			return new BTMutualInformationParameterSet(clazz, ess);
		}
	}

	/**
	 * Class for the parameters of a {@link BTMutualInformation} structure
	 * {@link Measure}.
	 * 
	 * @author Jan Grau
	 */
	public static class BTMutualInformationParameterSet extends
			InstanceParameterSet {

		/**
		 * Creates a new {@link BTMutualInformationParameterSet} with empty
		 * parameter values.
		 */
		public BTMutualInformationParameterSet() {
			super(BTMutualInformation.class);
		}

		/**
		 * Creates a new {@link BTMutualInformationParameterSet} with the
		 * parameter for the {@link DataSource} set to <code>clazz</code> and
		 * the parameter for the equivalent sample sizes (ess) set to
		 * <code>ess</code>.
		 * 
		 * @param clazz
		 *            the source of the data to compute the mutual information
		 * @param ess
		 *            the equivalent sample sizes for the foreground class and
		 *            the background, i.e. the background class or (in case of
		 *            more than two classes) all non-foreground classes
		 * 
		 * @throws Exception
		 *             if the parameters could not be created or set
		 */
		public BTMutualInformationParameterSet(DataSource clazz, double[] ess)
				throws Exception {
			super(BTMutualInformation.class);
			loadParameters();
			parameters.get(0).setValue(ess[0]);
			parameters.get(1).setValue(ess[1]);
			parameters.get(2).setValue(clazz);
		}

		/**
		 * The standard constructor for the interface {@link de.jstacs.Storable}
		 * . Creates a new {@link BTMutualInformationParameterSet} from its XML
		 * representation as defined by the {@link de.jstacs.Storable}
		 * interface.
		 * 
		 * @param representation
		 *            the XML code as {@link StringBuffer}
		 * 
		 * @throws NonParsableException
		 *             if the XML representation could not be parsed
		 */
		public BTMutualInformationParameterSet(StringBuffer representation)
				throws NonParsableException {
			super(representation);
		}

		/**
		 * Returns the equivalent sample sizes (ess) defined by this set of
		 * parameters.
		 * 
		 * @return the ess
		 */
		public double[] getEss() {
			return new double[] { (Double) parameters.get(0).getValue(),
					(Double) parameters.get(1).getValue() };
		}

		/**
		 * Returns the source of the data to compute the mutual information as
		 * defined by this set of parameters.
		 * 
		 * @return the source of the data
		 */
		public DataSource getClazz() {
			return (DataSource) ((EnumParameter) parameters.get(2)).getValue();
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see de.jstacs.parameters.InstanceParameterSet#getInstanceComment()
		 */
		@Override
		public String getInstanceComment() {
			return "Bayesian tree - mutual information";
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see de.jstacs.parameters.InstanceParameterSet#getInstanceName()
		 */
		@Override
		public String getInstanceName() {
			return "Build a Bayesian tree using mutual information as structure measure.";
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
			this.parameters.add(new EnumParameter(DataSource.class,
					"The data used to compute mutual information.", true));
		}

	}

}

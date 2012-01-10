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

package de.jstacs.differentiableStatisticalModels.directedGraphicalModels.structureLearning.measures.pmmMeasures;

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.algorithms.graphs.DAG;
import de.jstacs.algorithms.graphs.tensor.SymmetricTensor;
import de.jstacs.algorithms.graphs.tensor.Tensor;
import de.jstacs.data.DataSet;
import de.jstacs.differentiableStatisticalModels.directedGraphicalModels.structureLearning.measures.Measure;
import de.jstacs.differentiableStatisticalModels.directedGraphicalModels.structureLearning.measures.btMeasures.BTMutualInformation;
import de.jstacs.differentiableStatisticalModels.directedGraphicalModels.structureLearning.measures.btMeasures.BTMutualInformation.DataSource;
import de.jstacs.parameters.EnumParameter;
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.SimpleParameter;

/**
 * Class for the network structure of a
 * {@link de.jstacs.differentiableStatisticalModels.directedGraphicalModels.BayesianNetworkDiffSM}
 * that is a permuted Markov model based on mutual information.
 * 
 * @author Jan Grau
 * 
 */
public class PMMMutualInformation extends Measure {

	/**
	 * Creates a new {@link PMMMutualInformation} of order <code>order</code>.
	 * 
	 * @param order   the order, may be <code>1</code> or <code>2</code>.
	 * @param clazz    the classes used for computation of mutual information, as
	 *            defined by {@link DataSource}
	 * @param ess   the equivalent sample sizes of both classes
	 * 
	 * @throws Exception if the order is not <code>1</code> or <code>2</code>
	 */
	public PMMMutualInformation(byte order, DataSource clazz, double[] ess)
			throws Exception {
		this( new PMMMutualInformationParameterSet( order, clazz, ess ) );

	}

	/**
	 * Creates a new {@link PMMMutualInformation} from the corresponding
	 * {@link de.jstacs.parameters.InstanceParameterSet} <code>parameters</code>.
	 * 
	 * @param parameters  the corresponding parameters
	 * 
	 * @throws Exception   if the order is not <code>1</code> or <code>2</code>
	 */
	public PMMMutualInformation(PMMMutualInformationParameterSet parameters)
			throws Exception {
		super( parameters );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Recreates a {@link PMMMutualInformation} from its XML representation as
	 * returned by {@link #toXML()}.
	 * 
	 * @param buf  the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException if the XML code could not be parsed
	 */
	public PMMMutualInformation(StringBuffer buf) throws NonParsableException {
		super( buf );
	}	
	
	public String getXMLTag() {
		return "pmmMutualInformation";
	}

	/* (non-Javadoc)
	 * @see de.jstacs.differentiableStatisticalModels.directedGraphicalModels.structureLearning.measures.Measure#getInstanceName()
	 */
	@Override
	public String getInstanceName() {
		PMMMutualInformationParameterSet ps = (PMMMutualInformationParameterSet) parameters;
		String str = "Permuted Markov model of order " + ps.getOrder() + " with mutual information of";
		if (ps.getClazz() == DataSource.FG) {
			return str + " foreground";
		} else if (ps.getClazz() == DataSource.BG) {
			return str + " background";
		} else {
			return str + " foreground and background";
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.differentiableStatisticalModels.directedGraphicalModels.structureLearning.measures.Measure#getParents(de.jstacs.data.Sample, de.jstacs.data.Sample, double[], double[], int)
	 */
	@Override
	public int[][] getParents(DataSet fg, DataSet bg, double[] weightsFg,
			double[] weightsBg, int length) throws Exception {
		DataSet data = null;
		double[] weights = null;
		double ess2 = 0;
		PMMMutualInformationParameterSet ps = (PMMMutualInformationParameterSet) parameters;
		DataSource clazz = ps.getClazz();
		double[] ess = ps.getEss();
		byte order = ps.getOrder();
			
		if (clazz == DataSource.FG) {
			data = fg;
			weights = weightsFg;
			ess2 = ess[0];
		} else if (clazz == DataSource.BG) {
			data = bg;
			weights = weightsBg;
			ess2 = ess[1];
		} else {
			data = DataSet.union(fg, bg);
			weights = union(new double[][] { weightsFg, weightsBg });
			ess2 = ess[0] + ess[1];
		}

		Tensor t = new SymmetricTensor(length, order);
		fillTensor(t, getMI(getStatistics(data, weights, length, ess2),
				sum(weights) + ess2));
		if (order == 2) {
			fillTensor(t, getMI(getStatisticsOrderTwo(data, weights, length,
					ess2), sum(weights) + ess2));
		}
		int[] o = DAG.computeMaximalHP(t);

		return toParents(o, order);
	}

	/**
	 * Class for the parameters of a {@link PMMMutualInformation} structure
	 * {@link Measure}.
	 * 
	 * @author Jan Grau
	 */
	public static class PMMMutualInformationParameterSet extends
			MeasureParameterSet {

		/**
		 * Creates a new {@link PMMMutualInformationParameterSet} with empty
		 * parameter values.
		 * @throws ParameterException 
		 */
		public PMMMutualInformationParameterSet() throws ParameterException {
			super(BTMutualInformation.class);
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
			this.parameters
					.add(new SimpleParameter(
							DataType.BYTE,
							"Order",
							"The order of the permuted Markov model. Only 1 or 2 allowed.",
							true));
		}

		/**
		 * Creates a new {@link PMMMutualInformationParameterSet} with the
		 * parameter for the order set to <code>order</code>, the parameter for
		 * the {@link DataSource} set to <code>clazz</code> and the parameter
		 * for the equivalent sample sizes (ess) set to <code>ess</code>.
		 * 
		 * @param order  the order
		 * @param clazz  the source of the data to compute the mutual information
		 * @param ess  the equivalent sample sizes for the foreground class and
		 *            the background, i.e. the background class or (in case of
		 *            more than two classes) all non-foreground classes
		 *            
		 * @throws Exception if the parameters could not be created or set
		 */
		public PMMMutualInformationParameterSet(byte order, DataSource clazz,
				double[] ess) throws Exception {
			this();
			parameters.get(0).setValue(ess[0]);
			parameters.get(1).setValue(ess[1]);
			parameters.get(2).setValue(clazz);
			parameters.get(3).setValue(order);
		}

		/**
		 * Creates a new {@link PMMMutualInformationParameterSet} from its
		 * XML representation as defined by the {@link de.jstacs.Storable}
		 * interface.
		 * 
		 * @param representation the XML code as {@link StringBuffer}
		 * 
		 * @throws NonParsableException if the XML representation could not be parsed
		 */
		public PMMMutualInformationParameterSet(StringBuffer representation)
				throws NonParsableException {
			super(representation);
		}

		/**
		 * Returns the equivalent sample sizes (ess) defined by this set of parameters.
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

		/**
		 * Returns the order defined by this set of parameters.
		 * 
		 * @return the order
		 */
		public byte getOrder() {
			return (Byte) parameters.get(3).getValue();
		}

		/* (non-Javadoc)
		 * @see de.jstacs.parameters.InstanceParameterSet#getInstanceComment()
		 */
		@Override
		public String getInstanceComment() {
			return "Permuted Markov model - mutual information";
		}

		/* (non-Javadoc)
		 * @see de.jstacs.parameters.InstanceParameterSet#getInstanceName()
		 */
		@Override
		public String getInstanceName() {
			return "Build a permuted Markov model using mutual information as structure measure.";
		}
	}
}

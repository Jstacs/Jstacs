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

package de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.structureLearning.measures.pmmMeasures;

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.algorithms.graphs.DAG;
import de.jstacs.algorithms.graphs.tensor.SymmetricTensor;
import de.jstacs.algorithms.graphs.tensor.Tensor;
import de.jstacs.data.DataSet;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.DatatypeNotValidException;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.structureLearning.measures.Measure;

/**
 * Class for the network structure of a
 * {@link de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.BayesianNetworkDiffSM}
 * that is a permuted Markov model based on the explaining away residual.
 * 
 * @author Jan Grau
 */
public class PMMExplainingAwayResidual extends Measure {
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Recreates a {@link PMMExplainingAwayResidual} from its XML representation
	 * as returned by {@link #toXML()}.
	 * 
	 * @param buf
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML code could not be parsed
	 */
	public PMMExplainingAwayResidual(StringBuffer buf) throws NonParsableException {
		super( buf );
	}

	/**
	 * Creates a new {@link PMMExplainingAwayResidual} of order
	 * <code>order</code>.
	 * 
	 * @param order
	 *            the order, may be <code>1</code> or <code>2</code>
	 * @param ess
	 *            the equivalent sample sizes of both classes
	 * 
	 * @throws Exception
	 *             if the order is not <code>1</code> or <code>2</code>
	 */
	public PMMExplainingAwayResidual(byte order, double[] ess) throws Exception {
		this( new PMMExplainingAwayResidualParameterSet( order, ess ) );
	}

	/**
	 * Creates a new {@link PMMExplainingAwayResidual} from the corresponding
	 * {@link de.jstacs.parameters.InstanceParameterSet} <code>parameters/code>.
	 * 
	 * @param parameters
	 *            the corresponding parameters
	 * 
	 * @throws Exception
	 *             if the order is not <code>1</code> or <code>2</code>
	 */
	public PMMExplainingAwayResidual( PMMExplainingAwayResidualParameterSet parameters) throws Exception {
		super( parameters );
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.structureLearning.
	 * measures.Measure#getInstanceName()
	 */
	@Override
	public String getInstanceName() {
		return "Permuted Markov model of order " + ((PMMExplainingAwayResidualParameterSet)parameters).getOrder()
				+ " with explaining away residual";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.structureLearning.
	 * measures.Measure#getParents(de.jstacs.data.Sample, de.jstacs.data.Sample,
	 * double[], double[], int)
	 */
	@Override
	public int[][] getParents(DataSet fg, DataSet bg, double[] weightsFg,
			double[] weightsBg, int length) throws Exception {
		byte order = ((PMMExplainingAwayResidualParameterSet)parameters).getOrder();
		double[] ess = ((PMMExplainingAwayResidualParameterSet)parameters).getEss();
		
		Tensor t = new SymmetricTensor(length, order);

		double nFg = sum(weightsFg) + ess[0];
		double nBg = sum(weightsBg) + ess[1];

		fillTensor(t, getEAR(getStatistics(fg, weightsFg, length, ess[0]),
				getStatistics(bg, weightsBg, length, ess[1]), nFg, nBg));
		if (order == 2) {
			fillTensor(t, getEAR(getStatisticsOrderTwo(fg, weightsFg, length,
					ess[0]), getStatisticsOrderTwo(bg, weightsBg, length,
					ess[1]), nFg, nBg));
		}
		int[] o = DAG.computeMaximalHP(t);

		return toParents(o, order);
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.structureLearning.measures.Measure#getXMLTag()
	 */
	@Override
	public String getXMLTag() {
		return "pmmExplainingAwayResidual";
	}

	/**
	 * Class for the parameters of a {@link PMMExplainingAwayResidual} structure
	 * {@link Measure}.
	 * 
	 * @author Jan Grau
	 */
	public static class PMMExplainingAwayResidualParameterSet extends MeasureParameterSet {

		/**
		 * Creates a new {@link PMMExplainingAwayResidualParameterSet} with
		 * empty parameter values.
		 * @throws DatatypeNotValidException 
		 */
		public PMMExplainingAwayResidualParameterSet() throws DatatypeNotValidException {
			super(PMMExplainingAwayResidual.class);
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
			this.parameters
					.add(new SimpleParameter(
							DataType.BYTE,
							"Order",
							"The order of the permuted Markov model. Only 1 or 2 allowed.",
							true,
							new NumberValidator<Byte>((byte)1,(byte)2)));
		}

		/**
		 * Creates a new {@link PMMExplainingAwayResidualParameterSet} with the
		 * parameter for the order set to <code>order</code> and the parameter
		 * for the equivalent sample sizes (ess) set to <code>ess</code>.
		 * 
		 * @param order
		 *            the order
		 * @param ess
		 *            the equivalent sample sizes for the foreground class and
		 *            the background, i.e. the background class or (in case of
		 *            more than two classes) all non-foreground classes
		 * 
		 * @throws Exception
		 *             if the parameters could not be created or set
		 */
		public PMMExplainingAwayResidualParameterSet(byte order, double[] ess)
				throws Exception {
			this();
			parameters.get(0).setValue(ess[0]);
			parameters.get(1).setValue(ess[1]);
			parameters.get(2).setValue(order);
		}

		/**
		 * Creates a new {@link PMMExplainingAwayResidualParameterSet} from its
		 * XML representation as defined by the {@link de.jstacs.Storable}
		 * interface.
		 * 
		 * @param representation
		 *            the XML code as {@link StringBuffer}
		 * 
		 * @throws NonParsableException
		 *             if the XML representation could not be parsed
		 */
		public PMMExplainingAwayResidualParameterSet(StringBuffer representation)
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
		 * Returns the order defined by this set of parameters.
		 * 
		 * @return the order
		 */
		public byte getOrder() {
			return (Byte) parameters.get(2).getValue();
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see de.jstacs.parameters.InstanceParameterSet#getInstanceComment()
		 */
		@Override
		public String getInstanceComment() {
			return "Permuted Markov model - explaining away residual";
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see de.jstacs.parameters.InstanceParameterSet#getInstanceName()
		 */
		@Override
		public String getInstanceName() {
			return "Build a permuted Markov model using explaining away residual as structure measure.";
		}

	}

}

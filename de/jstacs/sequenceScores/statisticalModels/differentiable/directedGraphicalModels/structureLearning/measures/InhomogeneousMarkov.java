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

package de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.structureLearning.measures;

import de.jstacs.DataType;
import de.jstacs.data.DataSet;
import de.jstacs.io.NonParsableException;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.DatatypeNotValidException;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;

/**
 * Class for a network structure of a
 * {@link de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.BayesianNetworkDiffSM}
 * that is an inhomogeneous Markov model. The order of the Markov model can be
 * defined by the user. A Markov model of order <code>0</code> is also known as
 * position weight matrix (PWM), a Markov model of order <code>1</code> is also
 * known as weight array matrix (WAM) model.
 * 
 * @author Jan Grau
 */
public class InhomogeneousMarkov extends Measure {
	/**
	 * Creates the structure of an inhomogeneous Markov model of order
	 * <code>order</code>.
	 * 
	 * @param order
	 *            the order
	 * @throws CloneNotSupportedException if the parameters could not be cloned
	 * @throws IllegalValueException if the order is not allowed
	 */
	public InhomogeneousMarkov(int order) throws IllegalValueException, CloneNotSupportedException {
		this( new InhomogeneousMarkovParameterSet( order ) );
	}

	/**
	 * Creates a new {@link InhomogeneousMarkov} from the corresponding
	 * {@link de.jstacs.parameters.InstanceParameterSet} <code>parameters</code>.
	 * 
	 * @param parameters
	 *            the corresponding parameters
	 * @throws CloneNotSupportedException if the parameters could not be cloned
	 */
	public InhomogeneousMarkov(InhomogeneousMarkovParameterSet parameters) throws CloneNotSupportedException {
		super( parameters );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Recreates an {@link InhomogeneousMarkov} structure from its XML
	 * representation as returned by {@link #toXML()}.
	 * 
	 * @param buf
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML code could not be parsed
	 */
	public InhomogeneousMarkov(StringBuffer buf) throws NonParsableException {
		super( buf );
	}

	/**
	 * Returns the order of the Markov model as defined in the constructor
	 * 
	 * @return the order
	 */
	public int getOrder() {
		return ((InhomogeneousMarkovParameterSet)parameters).getOrder();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.structureLearning.
	 * measures.Measure#clone()
	 */
	@Override
	public InhomogeneousMarkov clone() throws CloneNotSupportedException {
		return (InhomogeneousMarkov) super.clone();
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
		return "Inhomogeneous Markov model of order " + getOrder();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.structureLearning.
	 * measures.Measure#getParents(de.jstacs.data.DataSet, de.jstacs.data.DataSet,
	 * double[], double[], int)
	 */
	@Override
	public int[][] getParents(DataSet fg, DataSet bg, double[] weightsFg,
			double[] weightsBg, int length) throws Exception {
		int[][] parents = new int[length][];
		int order = getOrder();
		for (int i = 0; i < parents.length; i++) {
			parents[i] = new int[(order < i ? order : i) + 1];
			for (int j = i; j >= i - order && j >= 0; j--) {
				parents[i][parents[i].length - (i - j) - 1] = j;
			}
		}
		return parents;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.structureLearning.
	 * measures.Measure#isShiftable()
	 */
	@Override
	public boolean isShiftable() {
		return true;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.structureLearning.measures.Measure#getXMLTag()
	 */
	@Override
	public String getXMLTag() {
		return "inhomogeneousMarkov";
	}

	/**
	 * Class for an {@link de.jstacs.parameters.InstanceParameterSet} that defines the parameters of
	 * an {@link InhomogeneousMarkov} structure {@link Measure}.
	 * 
	 * @author Jan Grau
	 */
	public static class InhomogeneousMarkovParameterSet extends MeasureParameterSet {

		/**
		 * Creates a new {@link InhomogeneousMarkovParameterSet} with empty
		 * parameter values.
		 */
		public InhomogeneousMarkovParameterSet() {
			super(InhomogeneousMarkov.class);
			try{
			this.parameters.add(new SimpleParameter(DataType.INT, "Order",
					"The order of the Markov model.", true));
			}catch(DatatypeNotValidException doesnothappen){
				throw new RuntimeException( doesnothappen );
			}
		}

		/**
		 * Creates a new {@link InhomogeneousMarkovParameterSet} with the
		 * parameter for the order set to <code>order</code>.
		 * 
		 * @param order
		 *            the order
		 * @throws IllegalValueException if the value of the order is not allowed
		 * 
		 */
		public InhomogeneousMarkovParameterSet(int order) throws IllegalValueException {
			this();
			this.parameters.get(0).setValue(order);
		}

		/**
		 * Creates a new {@link InhomogeneousMarkovParameterSet} from its XML
		 * representation as defined by the {@link de.jstacs.Storable}
		 * interface.
		 * 
		 * @param representation
		 *            the XML code as {@link StringBuffer}
		 * 
		 * @throws NonParsableException
		 *             if the XML representation could not be parsed
		 */
		public InhomogeneousMarkovParameterSet(StringBuffer representation)
				throws NonParsableException {
			super(representation);
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see de.jstacs.parameters.InstanceParameterSet#getInstanceComment()
		 */
		@Override
		public String getInstanceComment() {
			return "Inhomogeneous Markov model structure";
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see de.jstacs.parameters.InstanceParameterSet#getInstanceName()
		 */
		@Override
		public String getInstanceName() {
			return "Inhomogeneous Markov model";
		}

		/**
		 * Returns the order of the {@link InhomogeneousMarkov} structure
		 * measure as defined by this set of parameters.
		 * 
		 * @return the order
		 */
		public int getOrder() {
			return (Integer) parameters.get(0).getValue();
		}

	}

}

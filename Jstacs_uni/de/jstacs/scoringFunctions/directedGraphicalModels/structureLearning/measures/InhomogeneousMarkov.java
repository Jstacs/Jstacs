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

package de.jstacs.scoringFunctions.directedGraphicalModels.structureLearning.measures;

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.data.DataSet;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.InstanceParameterSet;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.DatatypeNotValidException;

/**
 * Class for a network structure of a
 * {@link de.jstacs.scoringFunctions.directedGraphicalModels.BayesianNetworkScoringFunction}
 * that is an inhomogeneous Markov model. The order of the Markov model can be
 * defined by the user. A Markov model of order <code>0</code> is also known as
 * position weight matrix (PWM), a Markov model of order <code>1</code> is also
 * known as weight array matrix (WAM) model.
 * 
 * @author Jan Grau
 */
public class InhomogeneousMarkov extends Measure {

	private int order;
	private InhomogeneousMarkovParameterSet parameters;

	/**
	 * Creates the structure of an inhomogeneous Markov model of order
	 * <code>order</code>.
	 * 
	 * @param order
	 *            the order
	 */
	public InhomogeneousMarkov(int order) {
		this.order = order;
	}

	/**
	 * Creates a new {@link InhomogeneousMarkov} from the corresponding
	 * {@link InstanceParameterSet} <code>parameters</code>.
	 * 
	 * @param parameters
	 *            the corresponding parameters
	 */
	public InhomogeneousMarkov(InhomogeneousMarkovParameterSet parameters) {
		this(parameters.getOrder());
		this.parameters = parameters;
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
		buf = XMLParser.extractForTag(buf, "inhomogeneousMarkov");
		this.order = XMLParser.extractObjectForTags(buf, "order", int.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
	}

	/**
	 * Returns the order of the Markov model as defined in the constructor
	 * 
	 * @return the order
	 */
	public int getOrder() {
		return order;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.scoringFunctions.directedGraphicalModels.structureLearning.
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
	 * de.jstacs.scoringFunctions.directedGraphicalModels.structureLearning.
	 * measures.Measure#getInstanceName()
	 */
	@Override
	public String getInstanceName() {
		return "Inhomogeneous Markov model of order " + order;
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
	public int[][] getParents(DataSet fg, DataSet bg, double[] weightsFg,
			double[] weightsBg, int length) throws Exception {
		int[][] parents = new int[length][];
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
	 * de.jstacs.scoringFunctions.directedGraphicalModels.structureLearning.
	 * measures.Measure#isShiftable()
	 */
	@Override
	public boolean isShiftable() {
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.Storable#toXML()
	 */
	public StringBuffer toXML() {
		StringBuffer buf = new StringBuffer();
		XMLParser.appendObjectWithTags(buf, order, "order");
		XMLParser.addTags(buf, "inhomogeneousMarkov");
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
			return new InhomogeneousMarkovParameterSet(order);
		}
	}

	/**
	 * Class for an {@link InstanceParameterSet} that defines the parameters of
	 * an {@link InhomogeneousMarkov} structure {@link Measure}.
	 * 
	 * @author Jan Grau
	 */
	public static class InhomogeneousMarkovParameterSet extends
			InstanceParameterSet {

		/**
		 * Creates a new {@link InhomogeneousMarkovParameterSet} with empty
		 * parameter values.
		 * @throws DatatypeNotValidException 
		 */
		public InhomogeneousMarkovParameterSet() throws DatatypeNotValidException {
			super(InhomogeneousMarkov.class);
			this.parameters.add(new SimpleParameter(DataType.INT, "Order",
					"The order of the Markov model.", true));
		}

		/**
		 * Creates a new {@link InhomogeneousMarkovParameterSet} with the
		 * parameter for the order set to <code>order</code>.
		 * 
		 * @param order
		 *            the order
		 * 
		 * @throws Exception
		 *             if the parameters could not be created or set
		 */
		public InhomogeneousMarkovParameterSet(int order) throws Exception {
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

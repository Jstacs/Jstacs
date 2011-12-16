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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.models.hmm.training;

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.algorithms.optimization.termination.AbstractTerminationCondition;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.validation.NumberValidator;


/**
 * This class implements an {@link de.jstacs.parameters.ParameterSet} for numerical training of an {@link de.jstacs.models.hmm.AbstractHMM}.
 * 
 * @author Jens Keilwagen
 */
public class NumericalHMMTrainingParameterSet extends MultiThreadedTrainingParameterSet {

	private static final String[] algorithmStrings = new String[]{	"steepest descent",
																	"conjugate gradients (F., R.)",
																	"conjugate gradients (P., R. positive)",
																	"quasi newton (D., F., P.)",
																	"quasi newton (B., F., G., S.)",
																	"limited memory quasi newton (B., F., G., S.; n=3)",
																	"limited memory quasi newton (B., F., G., S.; n=4)",
																	"limited memory quasi newton (B., F., G., S.; n=5)",
																	"limited memory quasi newton (B., F., G., S.; n=6)",
																	"limited memory quasi newton (B., F., G., S.; n=7)",
																	"limited memory quasi newton (B., F., G., S.; n=8)",
																	"limited memory quasi newton (B., F., G., S.; n=9)",
																	"limited memory quasi newton (B., F., G., S.; n=10)" };

	private static final Byte[] algorithms = new Byte[]{ Optimizer.STEEPEST_DESCENT,
														Optimizer.CONJUGATE_GRADIENTS_FR,
														Optimizer.CONJUGATE_GRADIENTS_PRP,
														Optimizer.QUASI_NEWTON_DFP,
														Optimizer.QUASI_NEWTON_BFGS,
														(byte)3,
														(byte)4,
														(byte)5,
														(byte)6,
														(byte)7,
														(byte)8,
														(byte)9,
														(byte)10 };
	
	/**
	 * This is the empty constructor that can be used to fill the parameters after creation.
	 */
	public NumericalHMMTrainingParameterSet() {
		addParameters();
	}

	/**
	 * This constructor can be used to create an instance with specified parameters.
	 * 
	 * @param starts the number of different starts
	 * @param tc the termination condition for stopping the algorithm
	 * @param threads the number of threads that should be used during optimization
	 * @param algorithm the algorithm that shall be used
	 * @param lineEps the threshold for stopping the line search
	 * @param startDist the start distance for the line search
	 * 
	 * @throws Exception if this {@link NumericalHMMTrainingParameterSet} could not be created
	 */
	public NumericalHMMTrainingParameterSet( int starts, AbstractTerminationCondition tc, int threads, byte algorithm, double lineEps, double startDist ) throws Exception {
		super( starts, tc, threads );
		addParameters();
		parameters.get( 3 ).setValue(algorithmStrings[getIndex( algorithmStrings, algorithms, algorithm, false )] );
		parameters.get( 4 ).setValue( lineEps );
		parameters.get( 5 ).setValue( startDist );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link NumericalHMMTrainingParameterSet} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link NumericalHMMTrainingParameterSet} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public NumericalHMMTrainingParameterSet( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}
	
	private void addParameters(){
		try {
			parameters.add( new SelectionParameter( DataType.BYTE,
					algorithmStrings,
					algorithms,
					"algorithm",
					"the algorithm that should be used for numerical optimization",
					true ) );

			parameters.add( new SimpleParameter( DataType.DOUBLE,
					"line epsilon",
					"the threshold for stopping the line search in the numerical training",
					true,
					new NumberValidator<Double>( 0d, Double.MAX_VALUE ) ) );
			parameters.add( new SimpleParameter( DataType.DOUBLE,
					"start distance",
					"the start distance for the line search in the numerical training",
					true,
					new NumberValidator<Double>( 0d, Double.MAX_VALUE ) ) );
		} catch ( ParameterException doesnothappen ) { } 
	}
	
	/**
	 * This method returns a byte encoding for the algorithm that should be used for optimization.
	 * 
	 * @return a byte encoding for the algorithm that should be used for optimization
	 */
	public byte getAlgorithm() {
		return (Byte) getParameterForName( "algorithm" ).getValue();
	}
	
	/**
	 * This method returns the threshold that should be used for stopping the line search during the optimization.
	 * 
	 * @return the threshold that should be used for stopping the line search during the optimization.
	 */
	public double getLineEps() {
		return (Double) getParameterForName( "line epsilon" ).getValue();
	}
	
	/**
	 * This method returns the start distance that should be used in the line search during the optimization.
	 * 
	 * @return the start distance that should be used in the line search during the optimization.
	 */
	public double getStartDistance() {
		return (Double) getParameterForName( "start distance" ).getValue();
	}
}

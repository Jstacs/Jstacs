/*
 * This file is part of Jstacs.
 * 
 * Jstacs is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 * 
 * Jstacs is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.sequenceScores.statisticalModels.trainable.discrete.inhomogeneous;

import java.io.IOException;
import java.io.OutputStream;
import java.util.Arrays;

import de.jstacs.algorithms.optimization.DifferentiableFunction;
import de.jstacs.algorithms.optimization.DimensionException;
import de.jstacs.algorithms.optimization.LimitedMedianStartDistance;
import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.algorithms.optimization.termination.TerminationCondition;
import de.jstacs.utils.RealTime;
import de.jstacs.utils.SafeOutputStream;
import de.jstacs.utils.Time;
import de.jstacs.utils.UserTime;

/**
 * 
 * @author Jens Keilwagen
 */
public class MEMTools {
	
	// algorithms for parameter estimation
	/**
	 * This constant can be used to specify that the model should use the iterative scaling for
	 * training.
	 */
	public static final byte SGIS_P = 11;

	/**
	 * This constant can be used to specify that the model should use the blockwise iterative scaling
	 * for training.
	 */
	public static final byte BGIS_P = 12;

	/**
	 * This constant can be used to specify that the model should use the iterative scaling for
	 * training.
	 */
	public static final byte GIS = 13;

	/**
	 * This constant can be used to specify that the model should use the iterative scaling for
	 * training.
	 */
	public static final byte SGIS = 14;

	/**
	 * This constant can be used to specify that the model should use the blockwise iterative scaling
	 * for training.
	 */
	public static final byte BGIS = 15;

	// others
	/**
	 * The epsilon for the line search in an optimization using the {@link Optimizer}.
	 */
	protected static final double LIN_EPS = 1E9;

	/**
	 * The start distance for the line search in an optimization using the {@link Optimizer}.
	 */
	protected static final double STARTDISTANCE = 1d;

	
	/**
	 * This method approximates the distribution either analytically or numerically.
	 * 
	 * @param constraints
	 *            the constraints to be used
	 * @param cond
	 *            the conditions
	 * @param sequence
	 *            the SequenceIterator used in normalization and numerical approximation
	 * @param algorithm
	 *            the choice of numerical approximation
	 * @param condition
	 *            the {@link TerminationCondition} for stopping the iterative algorithm
	 * @param stream
	 *            a possibility for writing some information
	 * @param alphLen
	 *            the alphabet length of all positions
	 * 
	 * @return the normalization constant
	 * 
	 * @throws Exception
	 *             if something went wrong inside the algorithms
	 */
	public static double train( MEMConstraint[] constraints, int[][] cond, SequenceIterator sequence, byte algorithm, TerminationCondition condition, OutputStream stream,
			int[] alphLen ) throws Exception
	{
		if( constraints.length == 0 ) {
			return 0;
		}
		SafeOutputStream sostream = SafeOutputStream.getSafeOutputStream(stream);
		setParametersToValue(constraints,0);

		double Z;
		int counter1 = 0;
		if( constraints.length > 1 )
		{
			for( int counter2 = 0; counter2 < constraints.length; counter2++ )
			{
				for( counter1 = 0; counter1 < constraints[counter2].getNumberOfSpecificConstraints(); counter1++ )
				{
					check( constraints[counter2].getFreq(counter1) );
				}
			}
			sostream.write( "numerically: " );
			Z = alphLen[0];
			for( counter1 = 1; counter1 < alphLen.length; counter1++ )
			{
				Z *= alphLen[counter1];
			}

			sequence.setBounds( alphLen );
			// chosen iterative algorithm
			switch( algorithm )
			{
			case SGIS_P:
				sostream.writeln( "SGIS_P" );
				counter1 = sequentialGeneralizedIterativeScalingP( constraints, Z, alphLen, sequence, condition, sostream );
				break;
			case BGIS_P:
				sostream.writeln( "BGIS_P" );
				counter1 = blockwiseGeneralizedIterativeScalingP( constraints, Z, alphLen, sequence, condition, sostream );
				break;
			case GIS:
				sostream.writeln( "GIS" );
				counter1 = generalizedIterativeScaling( constraints, Z, sequence, condition, sostream );
				break;
			case SGIS:
				sostream.writeln( "SGIS" );
				counter1 = sequentialGeneralizedIterativeScaling( constraints, Z, sequence, condition, sostream );
				break;
			case BGIS:
				sostream.writeln( "BGIS" );
				counter1 = blockwiseGeneralizedIterativeScaling( constraints, Z, sequence, condition, sostream );
				break;
			default:
				DualFunction psi = new DualFunction( sequence, constraints );
				double[] x = new double[psi.getDimensionOfScope()];
				counter1 = Optimizer.optimize( algorithm, psi, x, condition, LIN_EPS, new LimitedMedianStartDistance( 5, STARTDISTANCE ), sostream );
				psi.setValues( x );
				Z = psi.getZ();
			}
		}
		else
		{
			if( cond == null )
			{
				sostream.writeln( "unconditional analytically" );
				for( counter1 = 0; counter1 < constraints[0].getNumberOfSpecificConstraints(); counter1++ )
				{
					constraints[0].setExpLambda( counter1, constraints[0].getFreq( counter1 ) );
				}
			}
			else
			{
				sostream.writeln( "conditional analytically" );
				// createOffset and bins
				int[][] constrOffset = new int[cond.length][];
				double[][] constr = new double[cond.length][];
				sequence.setBounds( alphLen );
				int counter2, counter3, index;
				for( counter1 = 0; counter1 < cond.length; counter1++ )
				{
					constrOffset[counter1] = new int[cond[counter1].length];
					constrOffset[counter1][0] = 1;
					for( counter2 = 1; counter2 < cond[counter1].length; counter2++ )
					{
						constrOffset[counter1][counter2] = constrOffset[counter1][counter2 - 1] * alphLen[cond[counter1][counter2 - 1]];
					}
					constr[counter1] = new double[constrOffset[counter1][--counter2] * alphLen[cond[counter1][counter2]]];
				}
				// sum for subconstraints
				do
				{
					counter1 = constraints[0].satisfiesSpecificConstraint( sequence );
					for( counter2 = 0; counter2 < cond.length; counter2++ )
					{
						index = 0;
						for( counter3 = 0; counter3 < cond[counter2].length; counter3++ )
						{
							//index += constrOffset[counter2][counter3] * sequence.getValueAt( cond[counter2][counter3] );
							index += constrOffset[counter2][counter3] * sequence.seq[cond[counter2][counter3]];
						}
						constr[counter2][index] += constraints[0].getFreq( counter1 );
					}
				} while( sequence.next() );				
				// set values
				sequence.reset();
				double q;
				do
				{
					counter1 = constraints[0].satisfiesSpecificConstraint( sequence );
					q = constraints[0].getFreq( counter1 );
					for( counter2 = 0; counter2 < cond.length; counter2++ )
					{
						index = 0;
						for( counter3 = 0; counter3 < cond[counter2].length; counter3++ )
						{
							//index += constrOffset[counter2][counter3] * sequence.getValueAt( cond[counter2][counter3] );
							index += constrOffset[counter2][counter3] * sequence.seq[cond[counter2][counter3]];
						}
						q /= constr[counter2][index];
					}
					constraints[0].setExpLambda( counter1, q );
				} while( sequence.next() );
			}
			Z = 1;
			counter1 = 0;
		}
		return Z; 
	}

	/**
	 * This method is a convenience method that sets the same value for all parameter of the constraints
	 * @param constraint the constraints
	 * @param val the value to be set
	 */
	public static void setParametersToValue( MEMConstraint[] constraint, double val ) {
		for( int i = 0; i < constraint.length; i++ ) {
			for( int j = 0; j < constraint[i].getNumberOfSpecificConstraints(); j++ ) {
				constraint[i].setLambda(j, val);
			}
		}
	}
	
	private static double[][] createUniformP( int[] alphLen )
	{
		double anz = 1;
		int counter1 = 0, il1 = 1, il2 = 1;
		while( counter1 < alphLen.length && il2 < il2 * alphLen[counter1] )
		{
			anz *= alphLen[counter1];
			il2 *= alphLen[counter1++];
		}
		while( counter1 < alphLen.length && il1 < il1 * alphLen[counter1] )
		{
			anz *= alphLen[counter1];
			il1 *= alphLen[counter1++];
		}
		if( counter1 != alphLen.length )
		{
			throw new IllegalArgumentException( "Size of distribution too large." );
		}

		double[][] p = new double[il1][il2];
		double startvalue = 1d / anz;
		for( counter1 = 0; counter1 < il1; counter1++ )
		{
			Arrays.fill( p[counter1], startvalue );
		}
		return p;
	}

	/**
	 * This method computes the exponential part of the probability, i.e., everything except the normalization constant.
	 *  
	 * @param constraints the constraint
	 * @param fulfilled an array allowing to store which specific constraint is used, can be <code>null</code>
	 * @param sequence a sequence iteration
	 * 
	 * @return the exponential part of the probability
	 */
	public static double getExpPartOfProb( MEMConstraint[] constraints, int[] fulfilled, SequenceIterator sequence )
	{
		double p = 1;
		int idx;
		for( int counter = 0; counter < constraints.length; counter++ )
		{
			idx = constraints[counter].satisfiesSpecificConstraint( sequence );
			p *= constraints[counter].getExpLambda( idx );
			if( fulfilled != null ) {
				fulfilled[counter] = idx;
			}
		}
		return p;
	}

	private static void check( double p ) throws IllegalArgumentException
	{
		if( p <= 0 )
		{
			throw new IllegalArgumentException( "A marginal distribution became zero, please start the train-method again with a higher ESS." );
		}
	}
	
	private static Time getTimeObject() {
		Time t;
		try {
			t = new UserTime();
		} catch( Error e ) {
			System.out.println( "Warning: Could not load UserTime. Using RealTime instead." );
			t = new RealTime();
		}
		return t;
	}
	
	private static void out( SafeOutputStream sostream, int it, double time, double f, double delta ) throws IOException {
		sostream.writeln( it + " \t" + time + " \t" + f + " \t" + delta );
	}

	/**
	 * The blockwise generalized iterative scaling in the p-space by Jens Keilwagen.
	 */
	private static int blockwiseGeneralizedIterativeScalingP( MEMConstraint[] constraints, double Z, int[] alphLen, SequenceIterator sequence, TerminationCondition mode,
			SafeOutputStream sostream ) throws IOException, IllegalArgumentException
	{
		double[][] p = createUniformP(alphLen);
		int il1 = p.length, il2 = p[0].length;
		int counter1, counter2, counter3, iterations = 0, c = constraints.length - 1;
		double log_Z = Math.log( Z ), psi_old, psi_new = log_Z;
		double[] help_next, help_now = new double[constraints[0].getNumberOfSpecificConstraints()];

		Time t = getTimeObject();
		out(sostream,0,0,psi_new,0);
		// compute first marginal distribution
		sequence.reset();
		for( counter1 = 0; counter1 < il1; counter1++ )
		{
			for( counter2 = 0; counter2 < il2; counter2++ )
			{
				help_now[constraints[0].satisfiesSpecificConstraint( sequence )] += p[counter1][counter2];
				sequence.next();
			}
		}
		do
		{
			// main part
			for( counter3 = 0; counter3 < c; counter3++ )
			{
				// compute adjusting factor
				for( counter1 = 0; counter1 < help_now.length; counter1++ )
				{
					check( help_now[counter1] );
					help_now[counter1] = constraints[counter3].getFreq( counter1 ) / help_now[counter1];
					constraints[counter3].multiplyExpLambdaWith( counter1, help_now[counter1] );
				}
				// adjust
				sequence.reset();
				help_next = new double[constraints[counter3 + 1].getNumberOfSpecificConstraints()];
				for( counter1 = 0; counter1 < il1; counter1++ )
				{
					for( counter2 = 0; counter2 < il2; counter2++ )
					{
						p[counter1][counter2] *= help_now[constraints[counter3].satisfiesSpecificConstraint( sequence )];
						help_next[constraints[counter3 + 1].satisfiesSpecificConstraint( sequence )] += p[counter1][counter2];
						sequence.next();
					}
				}
				help_now = help_next;
			}
			// last adjustment
			for( counter1 = 0; counter1 < help_now.length; counter1++ )
			{
				check( help_now[counter1] );
				help_now[counter1] = constraints[c].getFreq( counter1 ) / help_now[counter1];
				constraints[c].multiplyExpLambdaWith( counter1, help_now[counter1] );
			}
			help_next = new double[constraints[0].getNumberOfSpecificConstraints()];
			sequence.reset();
			for( counter1 = 0; counter1 < il1; counter1++ )
			{
				for( counter2 = 0; counter2 < il2; counter2++ )
				{
					p[counter1][counter2] *= help_now[constraints[c].satisfiesSpecificConstraint( sequence )];
					help_next[constraints[0].satisfiesSpecificConstraint( sequence )] += p[counter1][counter2];
					sequence.next();
				}
			}
			help_now = help_next;

			psi_old = psi_new;
			psi_new = log_Z;
			for( counter1 = 0; counter1 < constraints.length; counter1++ )
			{
				for( counter2 = 0; counter2 < constraints[counter1].getNumberOfSpecificConstraints(); counter2++ )
				{
					psi_new -= constraints[counter1].getLambda( counter2 ) * constraints[counter1].getFreq( counter2 );
				}
			}
			out( sostream, ++iterations,  t.getElapsedTime(), psi_new, psi_old - psi_new );
		}
		while( mode.doNextIteration( iterations, psi_old, psi_new, null, null, 0, t ) );
		return iterations;
	}

	/**
	 * The blockwise generalized iterative scaling by Jens Keilwagen.
	 */
	private static int blockwiseGeneralizedIterativeScaling( MEMConstraint[] constraints, double Z, SequenceIterator sequence, TerminationCondition mode,
			SafeOutputStream sostream ) throws IOException, IllegalArgumentException
	{
		int counter1, counter2, iterations = 0, c = constraints.length - 1;
		double log_Z = Math.log( Z ), psi_old, psi_new = log_Z;
		double[] help_now = new double[constraints[0].getNumberOfSpecificConstraints()];

		Time t = getTimeObject();
		out(sostream,0,0,psi_new,0);
		do
		{
			// main part
			for( counter2 = 0; counter2 <= c; counter2++ )
			{
				sequence.reset();
				help_now = new double[constraints[counter2].getNumberOfSpecificConstraints()];
				do
				{
					help_now[constraints[counter2].satisfiesSpecificConstraint( sequence )] += getExpPartOfProb( constraints, null, sequence );
				}
				while( sequence.next() );
				// compute adjusting factor
				for( counter1 = 0; counter1 < help_now.length; counter1++ )
				{
					help_now[counter1] /= Z;
					check( help_now[counter1] );
					constraints[counter2].multiplyExpLambdaWith( counter1, constraints[counter2].getFreq( counter1 )
							/ help_now[counter1] );
				}
			}

			psi_old = psi_new;
			psi_new = log_Z;
			for( counter1 = 0; counter1 < constraints.length; counter1++ )
			{
				for( counter2 = 0; counter2 < constraints[counter1].getNumberOfSpecificConstraints(); counter2++ )
				{
					psi_new -= constraints[counter1].getLambda( counter2 ) * constraints[counter1].getFreq( counter2 );
				}
			}
			out( sostream, ++iterations,  t.getElapsedTime(), psi_new, psi_old - psi_new );
		}
		while( mode.doNextIteration( iterations, psi_old, psi_new, null, null, 0, t ) );
		return iterations;
	}

	/**
	 * The well know generalized iterative scaling algorithm by Darroch and Ratcliff.
	 */
	private static int generalizedIterativeScaling( MEMConstraint[] constraints, double Z, SequenceIterator sequence, TerminationCondition mode, SafeOutputStream sostream )
			throws IOException, IllegalArgumentException
	{
		int counter1, counter2, iterations = 0, n = constraints.length;
		double psi_new = Math.log( Z ), psi_old, p;
		double[][] adjust = new double[n][];
		for( counter1 = 0; counter1 < n; counter1++ )
		{
			adjust[counter1] = new double[constraints[counter1].getNumberOfSpecificConstraints()];
		}
		int[] fulfilled = new int[constraints.length];
		sequence.reset();
		Time t = getTimeObject();
		out(sostream,0,0,psi_new,0);
		do
		{
			p = getExpPartOfProb( constraints, fulfilled, sequence );
			for( counter1 = 0; counter1 < n; counter1++ )
			{
				adjust[counter1][fulfilled[counter1]] += p;
			}
		}
		while( sequence.next() );
		do
		{
			for( counter1 = 0; counter1 < n; counter1++ )
			{
				for( counter2 = 0; counter2 < adjust[counter1].length; counter2++ )
				{
					adjust[counter1][counter2] /= Z;
					check( adjust[counter1][counter2] );
					// change
					constraints[counter1].multiplyExpLambdaWith( counter2, Math.exp( Math.log( constraints[counter1]
							.getFreq( counter2 )
							/ adjust[counter1][counter2] )
							/ n ) );
					adjust[counter1][counter2] = 0;
				}
			}

			Z = 0;
			sequence.reset();
			do
			{
				p = getExpPartOfProb( constraints, fulfilled, sequence );
				Z += p;
				for( counter1 = 0; counter1 < n; counter1++ )
				{
					adjust[counter1][fulfilled[counter1]] += p;
				}
			}
			while( sequence.next() );

			psi_old = psi_new;
			psi_new = Math.log( Z );
			for( counter1 = 0; counter1 < constraints.length; counter1++ )
			{
				for( counter2 = 0; counter2 < constraints[counter1].getNumberOfSpecificConstraints(); counter2++ )
				{
					psi_new -= constraints[counter1].getLambda( counter2 ) * constraints[counter1].getFreq( counter2 );
				}
			}
			out( sostream, ++iterations,  t.getElapsedTime(), psi_new, psi_old - psi_new );
		}
		while( mode.doNextIteration( iterations, psi_old, psi_new, null, null, 0, t ) );
		return iterations;
	}

	/**
	 * The well known sequential generalized iterative scaling algorithm in the p-space.
	 */
	private static int sequentialGeneralizedIterativeScalingP( MEMConstraint[] constraints, double Z, int[] alphLen, SequenceIterator sequence, TerminationCondition mode,
			SafeOutputStream sostream ) throws IOException, IllegalArgumentException
	{
		double[][] p = createUniformP( alphLen );
		int il1 = p.length, il2 = p[0].length;
		int counter1, counter2, counter3, counter4, iterations = 0;
		double adjust_if, adjust_else, help;
		double psi_old, psi_new = Math.log( Z );

		Time t = getTimeObject();
		out(sostream,0,0,psi_new,0);
		do
		{
			// paper algorithm
			for( counter4 = 0; counter4 < constraints.length; counter4++ )
			{
				for( counter3 = 0; counter3 < constraints[counter4].getNumberOfSpecificConstraints(); counter3++ )
				{
					// compute marginal distribution
					help = 0;
					sequence.reset();
					for( counter1 = 0; counter1 < il1; counter1++ )
					{
						for( counter2 = 0; counter2 < il2; counter2++ )
						{
							if( constraints[counter4].satisfiesSpecificConstraint( sequence ) == counter3 )
							{
								help += p[counter1][counter2];
							}
							sequence.next();
						}
					}
					// adjusting
					check( help );
					adjust_if = constraints[counter4].getFreq( counter3 ) / help;
					adjust_else = (1 - constraints[counter4].getFreq( counter3 )) / (1 - help);
					constraints[counter4].multiplyExpLambdaWith( counter3, adjust_if * (1 - help)
							/ (1 - constraints[counter4].getFreq( counter3 )) );
					Z *= (1 - help) / (1 - constraints[counter4].getFreq( counter3 ));
					sequence.reset();
					for( counter1 = 0; counter1 < il1; counter1++ )
					{
						for( counter2 = 0; counter2 < il2; counter2++ )
						{
							if( constraints[counter4].satisfiesSpecificConstraint( sequence ) == counter3 )
							{
								p[counter1][counter2] *= adjust_if;
							}
							else
							{
								p[counter1][counter2] *= adjust_else;
							}
							sequence.next();
						}
					}
				}
			}
			// compute psi
			psi_old = psi_new;
			psi_new = Math.log( Z );
			for( counter1 = 0; counter1 < constraints.length; counter1++ )
			{
				for( counter2 = 0; counter2 < constraints[counter1].getNumberOfSpecificConstraints(); counter2++ )
				{
					psi_new -= constraints[counter1].getExpLambda( counter2 ) * constraints[counter1].getFreq( counter2 );
				}
			}
			out( sostream, ++iterations,  t.getElapsedTime(), psi_new, psi_old - psi_new );
		}
		while( mode.doNextIteration( iterations, psi_old, psi_new, null, null, 0, t ) );
		return iterations;
	}

	/**
	 * The well known sequential generalized iterative scaling algorithm.
	 */
	private static int sequentialGeneralizedIterativeScaling( MEMConstraint[] constraints, double Z, SequenceIterator sequence, TerminationCondition mode,
			SafeOutputStream sostream ) throws IOException, IllegalArgumentException
	{
		int counter1, counter2, counter3, iterations = 0, pos;
		double help, psi_old, psi_new = Math.log( Z );
		boolean isCorrect;

		Time t = getTimeObject();
		out(sostream,0,0,psi_new,0);
		do
		{
			// paper algorithm
			for( counter2 = 0; counter2 < constraints.length; counter2++ )
			{
				for( counter3 = 0; counter3 < constraints[counter2].getNumberOfSpecificConstraints(); counter3++ )
				{
					pos = constraints[counter2].getCorrectedPosition( 0 );
					// compute marginal distribution
					help = 0;
					sequence.reset();
					do
					{
						if( constraints[counter2].satisfiesSpecificConstraint( sequence ) == counter3 )
						{
							help += getExpPartOfProb( constraints, null, sequence );
							isCorrect = sequence.next();
						}
						else
						{
							isCorrect = sequence.skip( pos );
						}
					}
					while( isCorrect );
					// adjusting
					help /= Z;
					check( help );
					constraints[counter2].multiplyExpLambdaWith( counter3, constraints[counter2].getFreq( counter3 )
							* (1 - help) / (help * (1 - constraints[counter2].getFreq( counter3 ))) );
					Z *= (1 - help) / (1 - constraints[counter2].getFreq( counter3 ));
				}
			}
			// compute psi
			psi_old = psi_new;
			psi_new = Math.log( Z );
			for( counter1 = 0; counter1 < constraints.length; counter1++ )
			{
				for( counter2 = 0; counter2 < constraints[counter1].getNumberOfSpecificConstraints(); counter2++ )
				{
					psi_new -= constraints[counter1].getExpLambda( counter2 ) * constraints[counter1].getFreq( counter2 );
				}
			}
			out( sostream, ++iterations,  t.getElapsedTime(), psi_new, psi_old - psi_new );
		}
		while( mode.doNextIteration( iterations, psi_old, psi_new, null, null, 0, t ) );
		return iterations;
	}
	
	/**
	 * The dual function to the constraint problem of learning MEM's.
	 * 
	 * @author Jens Keilwagen
	 */
	public static class DualFunction extends DifferentiableFunction
	{
		private int n;

		private SequenceIterator it;

		private int[] shortcut;
		
		private MEMConstraint[] constraints;
		
		private double Z;

		/**
		 * The constructor of a dual function.
		 * 
		 * @param it
		 *            the correct initialized SequenceIterator
		 * @param constraints
		 *            the constraints used in this {@link DualFunction}
		 * 
		 * @throws IllegalArgumentException
		 *             if the constraints are not correct
		 */
		public DualFunction( SequenceIterator it, MEMConstraint[] constraints )
		{
			this.it = it;
			this.constraints = constraints;
			shortcut = new int[constraints.length];
			shortcut[0] = 0;
			n = constraints[0].getNumberOfSpecificConstraints();
			for( int i = 1; i < constraints.length; i++ )
			{
				shortcut[i] = shortcut[i - 1] + constraints[i - 1].getNumberOfSpecificConstraints();
				n += constraints[i].getNumberOfSpecificConstraints();
			}
		}

		public double evaluateFunction( double[] x ) throws DimensionException
		{
			if( x == null || x.length != getDimensionOfScope() )
			{
				if( x != null )
				{
					throw new DimensionException( x.length, n );
				}
				else
				{
					throw new DimensionException( 0, n );
				}
			}
			double[] y = exp( x );
			double erg = 0;
			it.reset();
			do
			{
				erg += getExpPartOfProb( y );
			}
			while( it.next() );
			erg = Math.log( erg );
			int counter1, counter2;
			for( counter1 = 0; counter1 < constraints.length; counter1++ )
			{
				for( counter2 = 0; counter2 < constraints[counter1].getNumberOfSpecificConstraints(); counter2++ )
				{
					erg -= constraints[counter1].getFreq( counter2 ) * x[shortcut[counter1] + counter2];
				}
			}
			return erg;
		}

		public double[] evaluateGradientOfFunction( double[] x ) throws DimensionException
		{
			double[] erg = new double[n], y = exp( x );
			double p, Z = 0;
			int counter1, counter2;
			int[] fulfilled = new int[constraints.length];

			it.reset();
			do
			{
				p = getExpPartOfProb( y, fulfilled );
				Z += p;
				for( counter1 = 0; counter1 < constraints.length; counter1++ )
				{
					erg[fulfilled[counter1]] += p;
				}
			}
			while( it.next() );

			for( counter1 = 0; counter1 < constraints.length; counter1++ )
			{
				for( counter2 = 0; counter2 < constraints[counter1].getNumberOfSpecificConstraints(); counter2++ )
				{
					erg[shortcut[counter1] + counter2] = erg[shortcut[counter1] + counter2] / Z
							- constraints[counter1].getFreq( counter2 );
				}
			}
			return erg;
		}

		public int getDimensionOfScope()
		{
			return n;
		}

		/**
		 * This method set the values of the Lagrange multiplicators of the constraints
		 * 
		 * @param x
		 *            the new values
		 */
		public void setValues( double[] x )
		{
			int counter1, counter2;
			double[] y = exp( x );
			for( counter1 = 0; counter1 < constraints.length; counter1++ )
			{
				for( counter2 = 0; counter2 < constraints[counter1].getNumberOfSpecificConstraints(); counter2++ )
				{
					constraints[counter1].setExpLambda( counter2, y[shortcut[counter1] + counter2] );
				}
			}

			Z = 0;
			it.reset();
			do
			{
				Z += getExpPartOfProb( y );
			}
			while( it.next() );
		}

		private double getExpPartOfProb( double[] x )
		{
			double erg = 1;
			int counter;
			for( counter = 0; counter < constraints.length; counter++ )
			{
				erg *= x[shortcut[counter] + constraints[counter].satisfiesSpecificConstraint( it )];
			}
			return erg;
		}

		private double getExpPartOfProb( double[] x, int[] fulfilled )
		{
			double erg = 1;
			int counter;
			for( counter = 0; counter < constraints.length; counter++ )
			{
				fulfilled[counter] = shortcut[counter] + constraints[counter].satisfiesSpecificConstraint( it );
				erg *= x[fulfilled[counter]];
			}
			return erg;
		}

		private double[] exp( double[] x )
		{
			double[] expX = new double[n];
			for( int counter = 0; counter < n; counter++ )
			{
				expX[counter] = Math.exp( x[counter] );
			}
			return expX;
		}
		
		private double getZ() {
			return Z;
		}
	}
}

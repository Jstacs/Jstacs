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

package de.jstacs.models.discrete.homogeneous;

import java.util.Random;

import de.jstacs.NonParsableException;
import de.jstacs.NotTrainedException;
import de.jstacs.WrongAlphabetException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.EmptySampleException;
import de.jstacs.data.Sample;
import de.jstacs.data.Sequence;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.WrongSequenceTypeException;
import de.jstacs.models.discrete.Constraint;
import de.jstacs.models.discrete.DGMParameterSet;
import de.jstacs.models.discrete.DiscreteGraphicalModel;
import de.jstacs.models.discrete.homogeneous.parameters.HomogeneousModelParameterSet;
import de.jstacs.results.NumericalResultSet;

/**
 * This class implements homogeneous models of arbitrary order.
 * 
 * @author Jens Keilwagen
 * 
 * @see HomogeneousModelParameterSet
 */
public abstract class HomogeneousModel extends DiscreteGraphicalModel {

	/**
	 * The powers of the alphabet length.
	 */
	protected int[] powers;

	/**
	 * The order of the model.
	 */
	protected byte order;

	/**
	 * Creates a homogeneous model from a parameter set.
	 * 
	 * @param params
	 *            the parameter set
	 * 
	 * @throws CloneNotSupportedException
	 *             if the parameter set could not be cloned
	 * @throws IllegalArgumentException
	 *             if the parameter set is not instantiated
	 * @throws NonParsableException
	 *             if the parameter set is not parsable
	 * 
	 * @see HomogeneousModelParameterSet
	 * @see DiscreteGraphicalModel#DiscreteGraphicalModel(DGMParameterSet)
	 */
	public HomogeneousModel( HomogeneousModelParameterSet params ) throws CloneNotSupportedException, IllegalArgumentException,
																	NonParsableException {
		super( params );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link HomogeneousModel} out of its XML representation.
	 * 
	 * @param stringBuff
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link HomogeneousModel} could not be reconstructed
	 *             out of the XML representation (the {@link StringBuffer} could
	 *             not be parsed)
	 * 
	 * @see de.jstacs.Storable
	 * @see DiscreteGraphicalModel#DiscreteGraphicalModel(StringBuffer)
	 */
	public HomogeneousModel( StringBuffer stringBuff ) throws NonParsableException {
		super( stringBuff );
	}

	/**
	 * Creates a {@link Sample} of a given number of {@link Sequence}s from a
	 * trained homogeneous model.
	 * 
	 * @param no
	 *            the number of {@link Sequence}s that should be in the
	 *            {@link Sample}
	 * @param length
	 *            the length of all {@link Sequence}s or an array of lengths
	 *            with the {@link Sequence} with index <code>i</code> having
	 *            length <code>length[i]</code>
	 * 
	 * @return the created {@link Sample}
	 * 
	 * @throws NotTrainedException
	 *             if the model was not trained
	 * @throws IllegalArgumentException
	 *             if the dimension of <code>length</code> is neither 1 nor
	 *             <code>no</code>
	 * @throws EmptySampleException
	 *             if <code>no == 0</code>
	 * @throws WrongSequenceTypeException
	 *             if the {@link Sequence} type is not suitable (for the
	 *             {@link AlphabetContainer})
	 * @throws WrongAlphabetException
	 *             if something is wrong with the alphabet
	 * 
	 * @see Sample#Sample(String, Sequence...)
	 */
	@Override
	public final Sample emitSample( int no, int... length ) throws NotTrainedException,
			IllegalArgumentException,
			EmptySampleException,
			WrongAlphabetException,
			WrongSequenceTypeException {
		if( !trained ) {
			throw new NotTrainedException();
		}
		Sequence[] seq = new Sequence[no];
		if( length.length == 1 ) {
			for( int i = 0; i < no; i++ ) {
				seq[i] = getRandomSequence( new Random(), length[0] );
			}
		} else if( length.length == no ) {
			for( int i = 0; i < no; i++ ) {
				seq[i] = getRandomSequence( new Random(), length[i] );
			}
		} else {
			throw new IllegalArgumentException( "The dimension of the array length is not correct." );
		}
		return new Sample( "sampled from " + getInstanceName(), seq );
	}

	/**
	 * This method creates a random {@link Sequence} from a trained homogeneous
	 * model.
	 * 
	 * @param r
	 *            the random generator
	 * @param length
	 *            the length of the {@link Sequence}
	 * 
	 * @return the created {@link Sequence}
	 * 
	 * @throws WrongSequenceTypeException
	 *             if the {@link Sequence} type is not suitable (for the
	 *             {@link AlphabetContainer})
	 * @throws WrongAlphabetException
	 *             if something is wrong with the alphabet
	 */
	protected abstract Sequence getRandomSequence( Random r, int length ) throws WrongAlphabetException, WrongSequenceTypeException;

	/* (non-Javadoc)
	 * @see de.jstacs.models.AbstractModel#getMaximalMarkovOrder()
	 */
	@Override
	public byte getMaximalMarkovOrder() {
		return order;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.Model#getNumericalCharacteristics()
	 */
	public NumericalResultSet getNumericalCharacteristics() throws Exception {
		return null;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.AbstractModel#getLogProbFor(de.jstacs.data.Sequence, int, int)
	 */
	@Override
	public final double getLogProbFor( Sequence sequence, int startpos, int endpos ) throws NotTrainedException, Exception {
		check( sequence, startpos, endpos );
		return logProbFor( sequence, startpos, endpos );
	}

	/**
	 * Trains the homogeneous model on all given {@link Sample}s.
	 * 
	 * @param data
	 *            the given {@link Sample}s
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see HomogeneousModel#train(Sample[], double[][])
	 */
	public void train( Sample[] data ) throws Exception {
		train( data, new double[data.length][] );
	}

	/**
	 * Trains the homogeneous model using an array of weighted {@link Sample}s.
	 * The {@link Sequence} weights in <code>weights[i]</code> are for the
	 * {@link Sample} in <code>data[i]</code>.
	 * 
	 * @param data
	 *            the given {@link Sample}s
	 * @param weights
	 *            the weights
	 * 
	 * @throws Exception
	 *             if something went wrong, furthermore <code>data.length</code>
	 *             has to be <code>weights.length</code>
	 */
	public abstract void train( Sample[] data, double[][] weights ) throws Exception;

	/* (non-Javadoc)
	 * @see de.jstacs.models.discrete.DiscreteGraphicalModel#set(de.jstacs.models.discrete.DGMParameterSet, boolean)
	 */
	@Override
	protected void set( DGMParameterSet params, boolean trained ) throws CloneNotSupportedException, NonParsableException {
		super.set( params, trained );
		order = (Byte)params.getParameterAt( 2 ).getValue();
		powers = new int[Math.max( order + 1, 2 )];
		powers[0] = 1;
		powers[1] = (int)alphabets.getAlphabetLengthAt( 0 );
		for( int i = 1; i < powers.length; i++ ) {
			powers[i] = powers[1] * powers[i - 1];
		}
	}

	/**
	 * Checks some constraints, these are in general conditions on the
	 * {@link de.jstacs.data.AlphabetContainer} of a (sub){@link Sequence}
	 * between <code>startpos</code> und <code>endpos</code>.
	 * 
	 * @param sequence
	 *            the {@link Sequence}
	 * @param startpos
	 *            the start position within the {@link Sequence}
	 * @param endpos
	 *            the end position within the {@link Sequence}
	 * 
	 * @throws NotTrainedException
	 *             if the model is not trained
	 * @throws IllegalArgumentException
	 *             if some arguments are wrong
	 * 
	 * @see DiscreteGraphicalModel#check(Sequence, int, int)
	 */
	@Override
	protected void check( Sequence sequence, int startpos, int endpos ) throws NotTrainedException, IllegalArgumentException {
		super.check( sequence, startpos, endpos );
		if( endpos >= sequence.getLength() ) {
			throw new IllegalArgumentException( "This endposition is impossible. Try: endposistion < sequence.length" );
		}
	}

	/**
	 * Chooses a value in <code>[0,end-start]</code> according to the
	 * distribution encoded in the frequencies of <code>distr</code> between the
	 * indices <code>start</code> and <code>end</code>.
	 * 
	 * <br>
	 * <br>
	 * 
	 * The instance <code>distr</code> is not changed in the process.
	 * 
	 * @param distr
	 *            the distribution
	 * @param start
	 *            the start index
	 * @param end
	 *            the end index
	 * @param randNo
	 *            a random number in [0,1]
	 * 
	 * @return the chosen value
	 * 
	 * @see Constraint#getFreq(int)
	 */
	protected final int chooseFromDistr( Constraint distr, int start, int end, double randNo ) {
		int c = start;
		while( randNo > distr.getFreq( c ) && c <= end ) {
			randNo -= distr.getFreq( c++ );
		}
		return c - start;
	}

	/**
	 * This method computes the logarithm of the probability of the given
	 * {@link Sequence} in the given interval. The method is only used in
	 * {@link de.jstacs.models.Model#getLogProbFor(Sequence, int, int)} after
	 * the method {@link HomogeneousModel#check(Sequence, int, int)} has been
	 * invoked.
	 * 
	 * @param sequence
	 *            the {@link Sequence}
	 * @param startpos
	 *            the start position within the {@link Sequence}
	 * @param endpos
	 *            the end position within the {@link Sequence}
	 * 
	 * @return the logarithm of the probability for the given subsequence
	 * 
	 * @see HomogeneousModel#check(Sequence, int, int)
	 * @see de.jstacs.models.Model#getLogProbFor(Sequence, int, int)
	 */
	protected abstract double logProbFor( Sequence sequence, int startpos, int endpos );

	/**
	 * Clones the given array of conditional probabilities.
	 * 
	 * @param p
	 *            the original conditional probabilities
	 * 
	 * @return an array of clones
	 */
	protected HomCondProb[] cloneHomProb( HomCondProb[] p ) {
		HomCondProb[] condProb = new HomCondProb[p.length];
		for( int i = 0; i < condProb.length; i++ ) {
			condProb[i] = new HomCondProb( p[i] );
		}
		return condProb;
	}

	/**
	 * This class handles the (conditional) probabilities of a homogeneous model
	 * in a fast way.
	 * 
	 * @author Jens Keilwagen
	 */
	protected class HomCondProb extends Constraint {

		private double[] lnFreq;

		/**
		 * The main constructor. Creates a new {@link HomCondProb} instance and
		 * checks that each position is used maximally once. In all/most cases
		 * <code>pos</code> is <code>new int[]{0,1,2...}</code> and <code>n</code>
		 * is <code>Math.pow({@link de.jstacs.data.Alphabet#length()},pos.length)</code>.
		 * 
		 * @param pos
		 *            the used positions (will be cloned), have to be
		 *            non-negative
		 * @param n
		 *            the number of specific constraints
		 * 
		 * @see Constraint#Constraint(int[], int)
		 */
		public HomCondProb( int[] pos, int n ) {
			super( pos, n );
		}

		/**
		 * The standard constructor for the interface {@link de.jstacs.Storable}
		 * . Creates a new {@link HomCondProb} instance out of its XML
		 * representation.
		 * 
		 * @param xml
		 *            the XML representation as {@link StringBuffer}
		 * 
		 * @throws NonParsableException
		 *             if the {@link StringBuffer} could not be parsed
		 * 
		 * @see de.jstacs.Storable
		 * @see Constraint#Constraint(StringBuffer)
		 */
		public HomCondProb( StringBuffer xml ) throws NonParsableException {
			super( xml );
		}

		/**
		 * Creates a new {@link HomCondProb} instance from a given one. This
		 * constructor is used for cloning instances, since any instance is an
		 * inner instance of a {@link HomogeneousModel}.
		 * 
		 * @param old
		 *            the old instance to be cloned
		 * 
		 * @see HomogeneousModel.HomCondProb#HomogeneousModel.HomCondProb(int[],int) HomogeneousModel.HomCondProb#HomCondProb(int[],int) 
		 */
		public HomCondProb( HomCondProb old ) {
			this( old.usedPositions, old.freq.length );
			System.arraycopy( old.freq, 0, freq, 0, freq.length );
			if( old.lnFreq != null ) {
				lnFreq( 0, freq.length );
			}
		}

		/* (non-Javadoc)
		 * @see de.jstacs.models.discrete.Constraint#estimate(double)
		 */
		@Override
		public void estimate( double ess ) {
			double pc = ess / (double)getNumberOfSpecificConstraints();
			if( usedPositions.length == 1 ) {
				estimateUnConditional( 0, freq.length, pc, false );
			} else {
				// conditional
				for( int counter1 = 0; counter1 < freq.length; counter1 += powers[1] ) {
					estimateUnConditional( counter1, counter1 + powers[1], pc, false );
				}
			}
		}

		/**
		 * Returns the logarithmic frequency at a given position
		 * <code>index</code>.
		 * 
		 * @param index
		 *            the given index
		 * 
		 * @return the logarithmic frequency at <code>index</code>
		 */
		public double getLnFreq( int index ) {
			return lnFreq[index];
		}

		/* (non-Javadoc)
		 * @see de.jstacs.models.discrete.Constraint#satisfiesSpecificConstraint(de.jstacs.data.Sequence, int)
		 */
		@Override
		public int satisfiesSpecificConstraint( Sequence seq, int start ) {
			int erg = 0, counter = 0, p = usedPositions.length - 1;
			for( ; counter < usedPositions.length; counter++, p-- ) {
				erg += powers[p] * seq.discreteVal( start + usedPositions[counter] );
			}
			return erg;
		}

		/* (non-Javadoc)
		 * @see de.jstacs.models.discrete.Constraint#toString()
		 */
		@Override
		public String toString() {
			String erg = "";
			int i = 1, l = usedPositions.length - 1;
			if( l > 0 ) {
				erg += usedPositions[0];
				while( i < l ) {
					erg += ", " + usedPositions[i++];
				}
				erg += " -> ";
			}
			return erg + usedPositions[l];
		}

		/**
		 * Adds the given <code>weight</code> to the counts corresponding to the
		 * {@link Sequence} <code>seq</code> from <code>start</code> to the end
		 * of the {@link Sequence}.
		 * 
		 * @param seq
		 *            the {@link Sequence}
		 * @param weight
		 *            the given weight
		 * @param start
		 *            the first index within the {@link Sequence}
		 * @param prevIndex
		 *            the previous index used for adding a count
		 */
		public final void addAll( Sequence seq, double weight, int start, int prevIndex ) {
			int l = seq.getLength();
			while( start < l ) {
				prevIndex = ( prevIndex % powers[order] ) * powers[1] + seq.discreteVal( start++ );
				counts[prevIndex] += weight;
			}

		}

		/* (non-Javadoc)
		 * @see de.jstacs.models.discrete.Constraint#appendAdditionalInfo(java.lang.StringBuffer)
		 */
		@Override
		protected void appendAdditionalInfo( StringBuffer xml ) {}

		private static final String XML_TAG = "HomCondProb";

		/* (non-Javadoc)
		 * @see de.jstacs.models.discrete.Constraint#getXMLTag()
		 */
		@Override
		protected String getXMLTag() {
			return XML_TAG;
		}

		/* (non-Javadoc)
		 * @see de.jstacs.models.discrete.Constraint#estimateUnConditional(int, int, double, boolean)
		 */
		@Override
		protected void estimateUnConditional( int start, int end, double pc, boolean exceptionWhenNoData ) {
			super.estimateUnConditional( start, end, pc, exceptionWhenNoData );
			lnFreq( start, end );
		}

		/**
		 * Computes the logarithm for all frequencies.
		 */
		private void lnFreq( int start, int end ) {
			if( lnFreq == null ) {
				lnFreq = new double[freq.length];
			}
			for( int i = start; i < end; i++ ) {
				lnFreq[i] = Math.log( freq[i] );
			}
		}

		/* (non-Javadoc)
		 * @see de.jstacs.models.discrete.Constraint#extractAdditionalInfo(java.lang.StringBuffer)
		 */
		@Override
		protected void extractAdditionalInfo( StringBuffer xml ) throws NonParsableException {
			lnFreq( 0, freq.length );
		}

		public String getDescription( AlphabetContainer con, int i ) {
			String res = null, s;
			DiscreteAlphabet d;
			for( int j = 0; j < usedPositions.length; j++ ) {
				d = (DiscreteAlphabet) con.getAlphabetAt( usedPositions[j] );
				s = "X_" + usedPositions[j] + "=" + d.getSymbolAt( i / powers[usedPositions.length-1-j] );
				if( res == null ) {
					res = s;
				} else {
					res = s + ", " + res;
				}
				i = i % powers[usedPositions.length-1-j];
			}
			res = res.replaceFirst( ", ", " | " );
			return "P(" + res + ")"; 
		}
	}
}

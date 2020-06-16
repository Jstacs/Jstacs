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
package de.jstacs.data.sequences;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.WrongLengthException;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;

/**
 * This class is for multidimensional sequences that can be used, for instance, for phylogenetic footprinting.
 * 
 * @author Jens Keilwagen
 * 
 * @param <T> the type of this {@link MultiDimensionalSequence} 
 */
public abstract class MultiDimensionalSequence<T> extends Sequence<T> {

	/**
	 * The internally used sequences.
	 */
	protected Sequence[] content;
	
	/**
	 * This constructor creates an {@link MultiDimensionalSequence} from a set of individual {@link Sequence}s.
	 * 
	 * @param seqAn the annotations for the aligned sequences
	 * @param sequence the individual sequences that have been aligned
	 * 
	 * @throws WrongLengthException if the sequences have different lengths
	 * @throws WrongAlphabetException if the sequences have different {@link de.jstacs.data.AlphabetContainer}s 
	 */
	public MultiDimensionalSequence( SequenceAnnotation[] seqAn, Sequence... sequence ) throws WrongLengthException, WrongAlphabetException {
		super( sequence[0].getAlphabetContainer(), seqAn );
		int l = sequence[0].getLength();
		for( int s = 1; s < sequence.length; s++ ) {
			if( sequence[s].getLength() != l ) {
				throw new WrongLengthException( "Creating an multi-dimensional sequence, all sequence have to have the same length" );
			}
			if( !sequence[s].getAlphabetContainer().checkConsistency( alphabetCon ) ) {
				throw new WrongAlphabetException( "All sequences have to have the same alphabet." );
			}
		}
		content = sequence.clone();
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.data.Sequence#continuousVal(int)
	 */
	@Override
	public double continuousVal(int pos) {
		return content[0].continuousVal( pos );
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.data.Sequence#discreteVal(int)
	 */
	@Override
	public int discreteVal(int pos) {
		return content[0].discreteVal( pos );
	}
	
	/**
	 * Returns a new instance of a {@link MultiDimensionalSequence} with given {@link SequenceAnnotation}s and given {@link Sequence}s.
	 * 
	 * @param seqAn the sequence annotations
	 * @param seqs the sequence buidl9ing the new {@link MultiDimensionalSequence}
	 * 
	 * @return a new instance of a {@link MultiDimensionalSequence} with given {@link SequenceAnnotation}s and given {@link Sequence}s
	 * 
	 * @throws WrongLengthException if the sequences have different lengths
	 * @throws WrongAlphabetException if the sequences have different {@link de.jstacs.data.AlphabetContainer}s
	 *  
	 * @see #flatCloneWithoutAnnotation()
	 * @see #complement()
	 * @see #reverseComplement()
	 */
	protected abstract MultiDimensionalSequence<T> getInstance( SequenceAnnotation[] seqAn, Sequence... seqs ) throws WrongLengthException, WrongAlphabetException;

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.data.Sequence#flatCloneWithoutAnnotation()
	 */
	@Override
	protected MultiDimensionalSequence<T> flatCloneWithoutAnnotation() {
		try {
			return getInstance( null, content );
		} catch ( WrongLengthException doesNotHappen ) {
			throw new RuntimeException( doesNotHappen.getMessage() );
		} catch ( WrongAlphabetException doesNotHappen ) {
			throw new RuntimeException( doesNotHappen.getMessage() );
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.data.Sequence#getEmptyRepresentation()
	 */
	@Override
	protected Object getEmptyRepresentation() {
		StringBuffer[] rep = new StringBuffer[content.length];
		for( int i = 0; i < content.length; i++ ) {
			rep[i] = new StringBuffer();
		}
		return rep;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.data.Sequence#addToRepresentation(java.lang.Object, int, java.lang.String)
	 */
	@Override
	protected void addToRepresentation( Object representation, int pos, String delim ) {
		for( int i = 0; i < content.length; i++ ) {
			((StringBuffer[])representation)[i].append( alphabetCon.getSymbol( pos, content[i].continuousVal( pos ) ) + delim );
		}
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.data.Sequence#getStringRepresentation(java.lang.Object)
	 */
	@Override
	protected String getStringRepresentation( Object representation ) {
		StringBuffer res = new StringBuffer();
		for( int i = 0; i < content.length; i++ ) {
			res.append( ((StringBuffer[])representation)[i] );
			res.append( "\n" );
		}
		return res.toString();
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.data.Sequence#getLength()
	 */
	@Override
	public int getLength() {
		return content[0].getLength();
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.data.Sequence#isMultiDimensional()
	 */
	@Override
	public boolean isMultiDimensional() {
		return true;
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.data.Sequence#complement(int, int)
	 */
	@Override
	public MultiDimensionalSequence<T> complement( int start, int end ) throws OperationNotSupportedException {
		if( alphabetCon.isReverseComplementable() ) {
			Sequence[] compContent = new Sequence[content.length];
			for( int s = 0; s < content.length; s++ ) {
				compContent[s] = content[s].complement( start, end );
			}
			try {
				return getInstance( null, compContent );
			} catch ( WrongLengthException doesNotHappen ) {
				throw new RuntimeException( doesNotHappen.getMessage() );
			} catch ( WrongAlphabetException doesNotHappen ) {
				throw new RuntimeException( doesNotHappen.getMessage() );
			}
		} else {
			throw new OperationNotSupportedException( "The alphabet of sequence has to be complementable." );
		}
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.data.Sequence#reverseComplement(int, int)
	 */
	@Override
	public MultiDimensionalSequence<T> reverseComplement( int start, int end ) throws OperationNotSupportedException {
		if( rc != null && start == 0 && end == getLength() ) {
			return (MultiDimensionalSequence) rc;
		} else if( alphabetCon.isReverseComplementable() ) {
			MultiDimensionalSequence<T> revComp;
			try {
				Sequence[] revCompContent = new Sequence[content.length];
				for( int s = 0; s < content.length; s++ ) {
					revCompContent[s] = content[s].reverseComplement( start, end );
				}
				revComp = getInstance( null, revCompContent );
				if( start == 0 && end == getLength() ) {
					rc = revComp;
					((MultiDimensionalSequence)rc).rc = this;
				}
				return revComp;
			} catch ( Exception e ) {
				// the current sequence is defined on the correct alphabet => so no WrongAlphabetException can occur
				RuntimeException doesNotHappen = new RuntimeException( e.getMessage() );
				doesNotHappen.setStackTrace( e.getStackTrace() );
				throw doesNotHappen;
			}
		} else {
			throw new OperationNotSupportedException( "The alphabet of sequence has to be reverse-complementable." );
		}
	}
	
	@Override
	public MultiDimensionalSequence<T> reverse( int start, int end ) throws OperationNotSupportedException {
		try {
			Sequence[] revContent = new Sequence[content.length];
			for( int s = 0; s < content.length; s++ ) {
				revContent[s] = content[s].reverse( start, end );
			}
			return getInstance( null, revContent );
		} catch( Exception e ) {
			OperationNotSupportedException onse = new OperationNotSupportedException(e.getMessage());
			onse.setStackTrace(e.getStackTrace());
			throw onse;
		}
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.data.Sequence#hashCodeForPos(int)
	 */
	@Override
	protected int hashCodeForPos( int pos ) {
		int h = 0;
		for( int i = 0; i < content.length; i++ ) {
			h = 31 * h + content[i].hashCodeForPos( pos );
		}
		return h;
	}
	
	/**
	 * This method returns the {@link SequenceAnnotation}[] for each dimension of this multidimensional sequence.
	 * 
	 * @return the {@link SequenceAnnotation} array for each dimension of this multidimensional sequence
	 * 
	 * @see Sequence#getAnnotation()
	 */
	public SequenceAnnotation[][] getAnnotations() {
		SequenceAnnotation[][] res = new SequenceAnnotation[content.length][];
		for(int i = 0; i < content.length; i++) {
				res[i] = content[i].getAnnotation();            
		}            
		return res;        
	}
	
	/**
	 * This method returns the number of internal sequences.
	 * 
	 * @return the number of internal sequences
	 * 
	 * @see #content
	 */
	public int getNumberOfSequences() {
		return content.length;
	}
	
	/**
	 * This method returns the internal sequence with index <code>index</code>.
	 * 
	 * @param index the index of the internal sequence
	 * 
	 * @return the internal sequence with index <code>index</code>
	 */
	public Sequence getSequence( int index ) {
		return content[index];
	}
}

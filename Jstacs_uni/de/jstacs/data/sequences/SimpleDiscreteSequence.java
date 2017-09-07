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

import java.util.Arrays;
import java.util.Random;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;

/**
 * This is the main class for any discrete sequence.
 * 
 * @author Jens Keilwagen
 */
public abstract class SimpleDiscreteSequence extends Sequence<int[]> {

	static Random r = new Random();
	
	/**
	 * This method implements the algorithm of D. Kandel et al. Discrete Applied Mathematics (1996) 171-185
	 * and returns a k-mer preserving shuffled sequence.
	 * 
	 * @param original the original sequence
	 * @param k the value for the k-mers, k should be larger than (or equal to) 1
	 * 
	 * @return the shuffled sequence
	 * 
	 * @throws Exception if the shuffled sequence could not be created
	 */
	public static SimpleDiscreteSequence shuffle( SimpleDiscreteSequence original, int k /*, boolean randomRotation*/ ) throws Exception {
		int n = original.getLength();
		if( n < 4*(k+1) ) {
			return original;
		}
		
		int[] shuffle = new int[n], help = shuffle.clone();
		for( int i = 0; i < n; i++ ) {
			shuffle[i] = original.discreteVal(i);
			//System.out.print( ( i % 10 == 0 ) ? "." : " " );
		}
		//System.out.println();
				 
		int anz = 0;
		for( int i = 0; i < shuffle.length; i++ ) {
			int a = r.nextInt(n-4*(k+1));
			int b = a+1+r.nextInt(n-3*(k+1)-a);
			int c = b+1+r.nextInt(n-2*(k+1)-b);
			int d = c+1+r.nextInt(n-k+1-c);
			
			boolean matches;
			int j = 0;
			//System.out.println();
			while( j < k-1 && shuffle[a+j] == shuffle[c+j] ) {
				//System.out.println( j + "\t" + (a+j) + "\t" + shuffle[a+j] + "\t" + (c+j) + "\t" + shuffle[c+j] );
				j++;
			}
			matches = j == k-1;
			j = 0;
			while( matches && j < k-1 && shuffle[b+j] == shuffle[d+j] ) {
				//System.out.println( j + "\t" + (b+j) + "\t" + shuffle[b+j] + "\t" + (d+j) + "\t" + shuffle[d+j] );
				j++;
			}
			matches &= j == k-1;
			if( matches ) {
				anz++;
				//System.out.println( i + "\t" + a + "\t" + b + "\t" + c + "\t" + d + "\t" + anz );
				Arrays.fill(help, 0);
				System.arraycopy( shuffle, 0, help, 0, a );
				System.arraycopy( shuffle, c, help, a, d-c );
				System.arraycopy( shuffle, b, help, a+d-c, c-b );
				System.arraycopy( shuffle, a, help, a+d-b, b-a );
				System.arraycopy( shuffle, d, help, d, n-d );
				System.arraycopy( help, 0, shuffle, 0, n );
			}
		}
		
		Exception ex = null;
		try {
			return new IntSequence( original.getAlphabetContainer(), shuffle );
		} catch( WrongAlphabetException wae ) {
			ex=wae;
		}  catch( WrongSequenceTypeException wste ) {
			ex=wste;
		}
		//will not happen.
		throw new RuntimeException(ex.getMessage());
	}
	
	/**
	 * This constructor creates a new {@link SimpleDiscreteSequence} with the
	 * {@link AlphabetContainer} <code>container</code> and the annotation
	 * <code>annotation</code> but without the content. The content has to be
	 * set by the constructor of the subclass.
	 * 
	 * @param container
	 *            the {@link AlphabetContainer} of the sequence
	 * @param annotation
	 *            the annotation of the sequence
	 * 
	 * @throws WrongAlphabetException
	 *             if the {@link AlphabetContainer} is not discrete
	 * 
	 * @see de.jstacs.data.sequences.Sequence#Sequence(AlphabetContainer, SequenceAnnotation[])
	 */
	public SimpleDiscreteSequence( AlphabetContainer container, SequenceAnnotation[] annotation ) throws WrongAlphabetException {
		super( container, annotation );
		if( !container.isDiscrete() ) {
			throw new WrongAlphabetException( "The alphabet is not discrete." );
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.Sequence#continuousVal(int)
	 */
	@Override
	public final double continuousVal( int pos ) {
		return discreteVal( pos );
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.data.Sequence#isMultiDimensional()
	 */
	public boolean isMultiDimensional()  {
		return false;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.data.Sequence#getEmptyContainer()
	 */
	public int[] getEmptyContainer() {
		return new int[1];
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.data.Sequence#fillContainer(java.lang.Object, int)
	 */
	public void fillContainer( int[] container, int pos ) {
		container[0] = discreteVal( pos );
	}
	
	protected Object getEmptyRepresentation() {
		return new StringBuffer();
	}
	
	protected void addToRepresentation( Object representation, int pos, String delim ) {
		((StringBuffer)representation).append( alphabetCon.getSymbol( pos, discreteVal( pos ) ) + delim );
	}
	
	protected String getStringRepresentation( Object representation ) {
		return representation.toString();
	}
	
	protected int hashCodeForPos( int pos ) {
		return discreteVal( pos );
	}
	
	public int compareTo( int[] t1, int[] t2 ) {
		if( t1.length == t2.length ) {
			for( int i = 0; i < t1.length; i++ ) {
				if( t1[i] != t2[i] ) {
					return t1[i]-t2[i];
				}
			}
			return 0;
		} else {
			return t1.length - t2.length;
		}
	}
}

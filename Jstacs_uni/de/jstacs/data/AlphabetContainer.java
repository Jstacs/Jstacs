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

package de.jstacs.data;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;

import de.jstacs.InstantiableFromParameterSet;
import de.jstacs.NonParsableException;
import de.jstacs.Storable;
import de.jstacs.WrongAlphabetException;
import de.jstacs.data.Alphabet.AlphabetParameterSet;
import de.jstacs.data.AlphabetContainerParameterSet.AlphabetArrayParameterSet;
import de.jstacs.data.AlphabetContainerParameterSet.SectionDefinedAlphabetParameterSet;
import de.jstacs.data.alphabets.ComplementableDiscreteAlphabet;
import de.jstacs.data.alphabets.ContinuousAlphabet;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.alphabets.DoubleSymbolException;
import de.jstacs.io.ParameterSetParser;
import de.jstacs.io.XMLParser;
import de.jstacs.io.ParameterSetParser.NotInstantiableException;
import de.jstacs.parameters.InstanceParameterSet;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.utils.SubclassFinder;

/**
 * The container for {@link Alphabet}s used in a {@link Sequence},
 * {@link Sample}, {@link de.jstacs.models.AbstractModel} or ... . The container
 * enables the user to have a different {@link Alphabet} at each position or at
 * least not the same {@link Alphabet} at all positions. This is impossible if
 * you use only instances of {@link Alphabet}. The container maps the given
 * {@link Alphabet} objects to the positions.
 * 
 * <br>
 * <br>
 * 
 * {@link AlphabetContainer} is immutable.
 * 
 * @author Jens Keilwagen
 * 
 * @see Alphabet
 */
public class AlphabetContainer implements Storable, InstantiableFromParameterSet, Comparable<AlphabetContainer> {

	private static final String XML_TAG = "AlphabetContainer";

	/**
	 * This <code>enum</code> defines types of {@link AlphabetContainer}s. These
	 * can be comprised of discrete alphabets, continuous alphabets, or
	 * discrete and continuous alphabets at different positions.
	 * 
	 * @author Jan Grau
	 */
	public enum AlphabetContainerType {
		/**
		 * The all {@link Alphabet}s are discrete.
		 */
		DISCRETE,
		/**
		 * The all {@link Alphabet}s are continuous.
		 */
		CONTINUOUS,
		/**
		 * The {@link Alphabet}s may be either continuous or discrete.
		 */
		BOTH;

		static AlphabetContainerType determineType( Alphabet[] alphabets ) {
			boolean discrete = true;
			boolean continuous = true;
			for( int i = 0; i < alphabets.length; i++ ) {
				discrete &= ( alphabets[i] instanceof DiscreteAlphabet );
				continuous &= ( alphabets[i] instanceof ContinuousAlphabet );
				if( !( discrete || continuous ) ) {
					return BOTH;
				}
			}
			AlphabetContainerType type = null;
			if( discrete ) {
				type = DISCRETE;
			} else {
				type = AlphabetContainerType.CONTINUOUS;
			}
			return type;
		}
		
		/**
		 * This method returns a {@link LinkedList} of
		 * {@link InstanceParameterSet}s which can be used to create
		 * {@link Alphabet}s that can be used in a {@link AlphabetContainer} of
		 * the given {@link AlphabetContainerType}.
		 * 
		 * @return a {@link LinkedList} of {@link InstanceParameterSet}s
		 * 
		 *         * @throws InstantiationException if any
		 *         {@link InstanceParameterSet} has no nullary constructor; or
		 *         if the instantiation fails for some other reason
		 * @throws IllegalAccessException
		 *             if any {@link InstanceParameterSet} or its nullary
		 *             constructor is not accessible
		 * @throws ClassNotFoundException
		 *             if one of the classes is present in the file system or
		 *             jar but cannot be loaded by the class loader
		 * @throws IOException
		 *             if the classes are searched for in a jar file, but that
		 *             file could not be accessed or read
		 */
		public LinkedList<InstanceParameterSet> getInstanceParameterSets()
				throws ClassNotFoundException, IOException,
				InstantiationException, IllegalAccessException {
			LinkedList<InstanceParameterSet> list = new LinkedList<InstanceParameterSet>();
			if (this != CONTINUOUS) {
				list.addAll(SubclassFinder.getInstanceParameterSets(
						DiscreteAlphabet.class, "de.jstacs.data"));
			}
			if (this != DISCRETE) {
				list
						.add(new ContinuousAlphabet.ContinuousAlphabetParameterSet());
			}
			return list;
		}
	}

	/**
	 * This method creates a new {@link AlphabetContainer} that uses as less as
	 * possible {@link Alphabet}s to describe the container. So, if possible,
	 * {@link Alphabet}s will be reused.
	 * 
	 * @param abc
	 *            the {@link Alphabet}s
	 * @param assignment
	 *            the assignment of the {@link Alphabet}s to the positions
	 * 
	 * @return an {@link AlphabetContainer} that uses as less as possible
	 *         {@link Alphabet}s
	 * 
	 * @see AlphabetContainer#AlphabetContainer(Alphabet[], int[])
	 */
	public static AlphabetContainer getSimplifiedAlphabetContainer( Alphabet[] abc, int[] assignment ) {
		ArrayList<Alphabet> list = new ArrayList<Alphabet>( abc.length );
		Alphabet current;
		int[] assign = new int[assignment.length], abcAssign = new int[abc.length];
		Arrays.fill( abcAssign, -1 );
		for( int j, i = 0; i < assign.length; i++ ) {
			if( abcAssign[assignment[i]] < 0 ) {
				j = 0;
				current = abc[assignment[i]];
				while( j < list.size() && !list.get( j ).checkConsistency( current ) ) {
					j++;
				}
				abcAssign[i] = j;
				if( j == list.size() ) {
					list.add( current );
				}
			}
			assign[i] = abcAssign[assignment[i]];
		}

		return new AlphabetContainer( list.toArray( new Alphabet[0] ), assign );
	}

	/**
	 * This method may be used to construct a new {@link AlphabetContainer} by
	 * incorporating additional {@link Alphabet}s into an existing
	 * {@link AlphabetContainer}.
	 * 
	 * @param aC
	 *            the {@link AlphabetContainer} used as template for the
	 *            returned {@link AlphabetContainer}
	 * @param a
	 *            the {@link Alphabet} that should be inserted
	 * @param useNewAlphabet
	 *            an array to define, which {@link Alphabet}s are used for which
	 *            positions <br>
	 *            <ul>
	 *            <li>the length of this array defines the length of the
	 *            {@link Sequence}s the returned {@link AlphabetContainer} is
	 *            capable to handle
	 *            <li>for each position <code>false</code> forces to use the
	 *            {@link Alphabet} given by the given {@link AlphabetContainer}
	 *            <li>for each position <code>true</code> forces to use the
	 *            given {@link Alphabet}
	 *            </ul>
	 *            Incorporation of the additional {@link Alphabet} must be
	 *            understood as defining new additional positions for which the
	 *            given {@link Alphabet} should be used. For example: let the
	 *            given {@link AlphabetContainer} contain three {@link Alphabet}
	 *            s <code>A0,A1,A2</code> (<code>Ai</code> for position
	 *            <code>i</code>) and therefore have a possible length of three.
	 *            Calling this method using this {@link AlphabetContainer}, an
	 *            additional {@link Alphabet} <code>A3</code> and an assignment
	 *            array of <code>[false</code>, <code>true</code>,
	 *            <code>true</code>, <code>false</code>, <code>false]</code>
	 *            returns a new {@link AlphabetContainer} having a possible
	 *            length of five and using the following alphabets for those
	 *            positions: <code>A0,A3,A3,A1,A2</code>. If the given
	 *            {@link AlphabetContainer} has a possible length not equal
	 *            zero, then the assignment array must contain as many
	 *            <code>false</code>-values as the length of the given
	 *            {@link AlphabetContainer}.
	 * 
	 * @return a new {@link AlphabetContainer} as described above
	 * 
	 * @throws IllegalArgumentException
	 *             if <code>useNewAlphabet</code> is <code>null</code> or has
	 *             length 0
	 * 
	 * @see AlphabetContainer#AlphabetContainer(Alphabet[], int[])
	 */
	public static AlphabetContainer insertAlphabet( AlphabetContainer aC, Alphabet a, boolean[] useNewAlphabet ) throws IllegalArgumentException {

		if( useNewAlphabet == null || useNewAlphabet.length == 0 ) {
			throw new IllegalArgumentException( "given useNewAlphabet-array is null or has length 0" );
		}

		Alphabet[] tempAs = new Alphabet[aC.alphabet.length + 1];
		tempAs[tempAs.length - 1] = a;
		int[] tempAssign = new int[useNewAlphabet.length];

		int i;
		int pos1 = 0;

		if( aC.getPossibleLength() == 0 ) {

			for( i = 0; i < tempAssign.length; i++ ) {

				if( useNewAlphabet[i] ) {
					tempAssign[i] = 0;
				} else {
					tempAssign[i] = 1;
				}
			}

		} else {

			for( i = 0; i < tempAssign.length; i++ ) {

				if( useNewAlphabet[i] ) {
					tempAssign[i] = aC.index[pos1++];
				} else {
					tempAssign[i] = tempAs.length - 1;
				}
			}
		}
		return new AlphabetContainer( tempAs, tempAssign );
	}

	/**
	 * The underlying {@link Alphabet}s of this {@link AlphabetContainer}.
	 */
	private Alphabet[] alphabet;

	private String delim;

	/**
	 * The array with the assignment of the {@link Alphabet}s to the positions
	 * (in the inhomogeneous case).
	 */
	private int[] index;

	private AlphabetContainerParameterSet parameters;

	private double l;

	/**
	 * Creates a new simple {@link AlphabetContainer}. All positions use the
	 * same {@link Alphabet} and therefore sequences of arbitrary length can be
	 * handled.
	 * 
	 * @param abc
	 *            the {@link Alphabet} for all positions
	 */
	public AlphabetContainer( Alphabet abc ) {
		index = null;
		alphabet = new Alphabet[]{ abc };
		precompute();
	}

	/**
	 * Creates a new {@link AlphabetContainer} with different {@link Alphabet}s
	 * for each position. The assignment of the {@link Alphabet}s to the
	 * positions is given by the order in the {@link Alphabet} array. This
	 * constructor should only be used if all {@link Alphabet}s are pairwise
	 * different.
	 * 
	 * @param abc
	 *            the different {@link Alphabet}s for each position
	 * 
	 * @see AlphabetContainer#AlphabetContainer(Alphabet[], int[])
	 */
	public AlphabetContainer( Alphabet... abc ) {
		this( abc, null );
	}

	/**
	 * Creates an new sparse {@link AlphabetContainer} based on given
	 * {@link AlphabetContainer}s.
	 * 
	 * @param cons
	 *            the given {@link AlphabetContainer}s
	 * @param lengths
	 *            the corresponding lengths of each {@link AlphabetContainer}
	 *            that is used
	 * 
	 * @throws IllegalArgumentException
	 *             if the given length for an {@link AlphabetContainer} is not
	 *             possible
	 */
	public AlphabetContainer( AlphabetContainer[] cons, int[] lengths ) throws IllegalArgumentException {
		int j, i = 0, n = 0, k;
		for( ; i < lengths.length; i++ ) {
			j = cons[i].getPossibleLength();
			if( j != 0 && j != lengths[i] ) {
				throw new IllegalArgumentException( "The AlphabetContainer " + i
													+ " is not able to handle sequences of length "
													+ lengths[i]
													+ "." );
			} else {
				n += lengths[i];
			}
		}
		index = new int[n];
		int[] help;
		ArrayList<Alphabet> abcs = new ArrayList<Alphabet>( n );
		for( n = i = 0; i < lengths.length; i++ ) {
			help = new int[cons[i].alphabet.length];
			for( j = 0; j < help.length; j++ ) {
				k = 0;
				while( k < abcs.size() && !cons[i].alphabet[j].checkConsistency( abcs.get( k ) ) ) {
					k++;
				}
				if( k == abcs.size() ) {
					// no consistent alphabet was found
					// -> the alphabet need to be added
					abcs.add( cons[i].alphabet[j] );
				}
				// else // a consistent alphabet was found
				help[j] = k;
			}
			if( cons[i].index == null ) {
				for( j = 0; j < lengths[i]; j++ ) {
					index[n++] = help[0];
				}
			} else {
				for( j = 0; j < lengths[i]; j++ ) {
					index[n++] = help[cons[i].index[j]];
				}
			}
		}
		alphabet = abcs.toArray( new Alphabet[0] );
		if( alphabet.length == 1 ) {
			index = null;
		}
	}

	/**
	 * Creates a new {@link AlphabetContainer} that uses different
	 * {@link Alphabet}s. The {@link Alphabet}s can be used more than once. The
	 * assignment for the {@link Alphabet}s to the positions is given by the
	 * array <code>assignment</code>.
	 * 
	 * @param abc
	 *            the {@link Alphabet}s
	 * @param assignment
	 *            the assignment array
	 * 
	 * @throws IllegalArgumentException
	 *             if the assignment of the {@link Alphabet}s to the positions
	 *             is not correct
	 */
	public AlphabetContainer( Alphabet[] abc, int[] assignment ) throws IllegalArgumentException {
		if( abc.length == 1 ) {
			index = null;
		} else {
			int i = 0;
			if( assignment == null ) {
				index = new int[abc.length];
				while( i < index.length ) {
					index[i] = i++;
				}
			} else {
				index = new int[assignment.length];
				boolean[] used = new boolean[abc.length];
				while( i < index.length && 0 <= assignment[i] && assignment[i] < abc.length ) {
					used[assignment[i]] = true;
					index[i] = assignment[i++];
				}
				if( i < index.length ) {
					throw new IllegalArgumentException( "The assignment from the positions to the alphabets is corrupted at position " + i
														+ "." );
				}
				i = 1;
				while( i < used.length && ( used[0] &= used[i++] ) );
				if( !used[0] ) {
					throw new IllegalArgumentException( "The assignment from the positions to the alphabets is not surjective. (Not all alphabets are used.)" );
				}
			}
		}
		alphabet = abc.clone();
		precompute();
	}

	/**
	 * Creates a new {@link AlphabetContainer} from an
	 * {@link AlphabetContainerParameterSet} that contains all necessary
	 * parameters.
	 * 
	 * @param parameters
	 *            the parameter set
	 * 
	 * @throws IllegalArgumentException
	 *             if something is wrong with the parameters in the
	 *             {@link AlphabetContainerParameterSet}
	 * @throws DoubleSymbolException
	 *             if the definitions within <code>parameters</code> contains a
	 *             symbol twice
	 * @throws NotInstantiableException
	 *             if an instance could not be created
	 */
	public AlphabetContainer( AlphabetContainerParameterSet parameters ) throws IllegalArgumentException, DoubleSymbolException,
																		NotInstantiableException {
		try {
			this.parameters = parameters.clone();
			ParameterSet alphSet = (ParameterSet)parameters.getParameterAt( 0 ).getValue();
			if( alphSet instanceof AlphabetParameterSet ) {
				// index = null;
				index = new int[1];
				alphabet = new Alphabet[]{ (Alphabet)ParameterSetParser.getInstanceFromParameterSet( (InstanceParameterSet)alphSet ) };
			} else if( alphSet instanceof AlphabetArrayParameterSet ) {
				index = new int[(Integer)alphSet.getParameterAt( 0 ).getValue()];
				alphabet = new Alphabet[alphSet.getNumberOfParameters() - 1];
				for( int i = 0; i < index.length; i++ ) {
					index[i] = i;
					ParameterSet par = (ParameterSet)alphSet.getParameterAt( i + 1 ).getValue();
					if( par instanceof SimpleParameterSet ) {
						par = (ParameterSet)par.getParameterAt( 0 ).getValue();
					}
					if( par instanceof AlphabetParameterSet ) {
						alphabet[i] = (Alphabet)ParameterSetParser.getInstanceFromParameterSet( (InstanceParameterSet)par );
					} else {
						throw new IllegalArgumentException( "The parameters of the alphabet at " + i + " are not given correctly" );
					}
				}
			} else if( alphSet instanceof SectionDefinedAlphabetParameterSet ) {
				if( !alphSet.hasDefaultOrIsSet() ) {
					throw new IllegalArgumentException( alphSet.getErrorMessage() );
				}
				index = new int[(Integer)alphSet.getParameterAt( 0 ).getValue()];
				alphabet = new Alphabet[alphSet.getNumberOfParameters() - 1];
				for( int i = 0; i < alphabet.length; i++ ) {
					ParameterSet singleAlph = (ParameterSet)alphSet.getParameterAt( i + 1 ).getValue();

					if( singleAlph.getParameterAt( 0 ).getValue() instanceof AlphabetParameterSet ) {
						alphabet[i] = (Alphabet)ParameterSetParser.getInstanceFromParameterSet( (InstanceParameterSet)singleAlph.getParameterAt( 0 )
								.getValue() );
					} else {
						throw new IllegalArgumentException( "Parameter for alphabet no. " + i
															+ " of unexpected type: "
															+ singleAlph.getClass()
															+ "." );
					}
					String section = (String)singleAlph.getParameterAt( 1 ).getValue();
					try {
						Iterator<Integer> posIt = SectionDefinedAlphabetParameterSet.parseSections( section ).iterator();
						while( posIt.hasNext() ) {
							index[posIt.next()] = i;
						}
					} catch ( Exception e ) {
						e.printStackTrace();
						throw new IllegalArgumentException( "Malformed sections for alphabet no. " + i + "." );
					}
				}
			} else {
				throw new IllegalArgumentException( "Wrong parameter type" );
			}
			precompute();
		} catch ( CloneNotSupportedException e ) {
			throw new IllegalArgumentException( e.getCause().getMessage() );
		}
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link AlphabetContainer} out of its XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link AlphabetContainer} could not be reconstructed
	 *             out of the XML representation (the {@link StringBuffer} could
	 *             not be parsed)
	 * 
	 * @see de.jstacs.Storable
	 */
	public AlphabetContainer( StringBuffer xml ) throws NonParsableException {
		StringBuffer buf = XMLParser.extractForTag( xml, XML_TAG );
		alphabet = XMLParser.extractObjectForTags( buf, "Alphabets", Alphabet[].class );

		if( alphabet.length > 1 ) {
			index = XMLParser.extractObjectForTags( buf, "Assignment", int[].class );
			int i = 0;
			while( i < index.length && index[i] < alphabet.length ) {
				i++;
			}
			if( i < index.length ) {
				throw new NonParsableException( "The assignment from the positions to the alphabets is corrupted at position " + i + "." );
			}
		}
		precompute();
	}

	/**
	 * Checks if this {@link AlphabetContainer} is consistent consistent with
	 * another {@link AlphabetContainer}.
	 * 
	 * @param abc
	 *            the second {@link AlphabetContainer}
	 * 
	 * @return <code>true</code> if the {@link AlphabetContainer}s are
	 *         consistent, <code>false</code> otherwise
	 * 
	 * @see Comparable#compareTo(Object)
	 */
	public boolean checkConsistency( AlphabetContainer abc ) {
		return compareTo( abc ) == 0;
	}

	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	public int compareTo( AlphabetContainer abc ) {
		if( abc == this ) {
			return 0;
		}
		int i = 0, a1 = getPossibleLength(), a2 = abc.getPossibleLength();
		if( a1 == a2 ) {
			if( a1 == 0 ) {
				return alphabet[0].compareTo( abc.alphabet[0] );
			} else {
				boolean erg = true;
				// to be quicker
				byte notChecked = 0, consistent = 1, inconsistent = 2;
				byte[][] checked = new byte[alphabet.length][abc.alphabet.length];
				int current = 0;
				while( i < a2 && erg ) {
					if( checked[index[i]][abc.index[i]] == notChecked ) {
						current = alphabet[index[i]].compareTo( abc.alphabet[abc.index[i]] );
						checked[index[i]][abc.index[i]] = current == 0 ? consistent : inconsistent;
						erg &= checked[index[i]][abc.index[i]] == consistent;
					}
					i++;
				}
				return current;
			}
		} else {
			return a1 - a2;
		}
	}

	private void precompute() {
		int i = 0;
		while( i < alphabet.length && alphabet[i] instanceof DiscreteAlphabet
				&& ( (DiscreteAlphabet)alphabet[i] ).getMaximalSymbolLength() == 1 ) {
			i++;
		}
		delim = ( i == alphabet.length ) ? "" : " ";
		l = alphabet[0].length();
		for( i = 1; i < alphabet.length; i++ ) {
			l = Math.max( l, alphabet[i].length() );
		}
	}

	/**
	 * Returns the underlying {@link Alphabet} of position <code>pos</code>.
	 * Please note that the {@link Alphabet} is returned as reference, so take
	 * care of what you are doing with it!
	 * 
	 * @param pos
	 *            the position
	 * 
	 * @return the {@link Alphabet} of the given position
	 */
	public Alphabet getAlphabetAt( int pos ) {
		return alphabet[getAlphabetIndexForPosition( pos )];
	}

	/**
	 * Returns the length of the underlying {@link Alphabet} of position
	 * <code>pos</code>.
	 * 
	 * @param pos
	 *            the position
	 * 
	 * @return the length of the underlying {@link Alphabet} of position
	 *         <code>pos</code>
	 */
	public double getAlphabetLengthAt( int pos ) {
		return getAlphabetAt( pos ).length();
	}

	/**
	 * Returns the encoded symbol for <code>sym</code> of the {@link Alphabet}
	 * of position <code>pos</code> of this {@link AlphabetContainer}.
	 * 
	 * @param pos
	 *            the position of the {@link Alphabet}
	 * @param sym
	 *            the symbol that should be returned encoded
	 * 
	 * @return the encoded symbol
	 * 
	 * @throws WrongAlphabetException
	 *             if the symbol is not defined in the {@link Alphabet} of the
	 *             given position
	 * 
	 * @see AlphabetContainer#getAlphabetAt(int)
	 */
	public double getCode( int pos, String sym ) throws WrongAlphabetException {
		if( isDiscreteAt( pos ) ) {
			return ( (DiscreteAlphabet)getAlphabetAt( pos ) ).getCode( sym );
		} else {
			double candidat = Double.parseDouble( sym );
			if( !( (ContinuousAlphabet)getAlphabetAt( pos ) ).isEncodedSymbol( candidat ) ) {
				throw new WrongAlphabetException();
			}
			return candidat;
		}
	}

	/**
	 * Returns an {@link AlphabetContainer} of {@link Alphabet}s e.g. for
	 * composite motifs/sequences.
	 * 
	 * @param start
	 *            the array of start indices
	 * @param length
	 *            the array of lengths
	 * 
	 * @return the {@link AlphabetContainer}
	 * 
	 * @see AlphabetContainer#getSubContainer(int, int)
	 */
	public AlphabetContainer getCompositeContainer( int[] start, int[] length ) {
		if( alphabet.length == 1 ) {
			return this;
		} else {
			int i = 0, j, l = 0;
			for( ; i < length.length; i++ ) {
				l += length[i];
			}
			int[] ind = new int[l];
			int[] used = new int[alphabet.length];
			Arrays.fill( used, -1 );
			ArrayList<Alphabet> list = new ArrayList<Alphabet>();
			l = 0;
			for( i = 0; i < length.length; i++ ) {
				for( j = 0; j < length[i]; j++ ) {
					if( used[index[start[i] + j]] < 0 ) {
						used[index[start[i] + j]] = list.size();
						list.add( alphabet[index[start[i] + j]] );
					}
					ind[l++] = used[index[start[i] + j]];
				}
			}
			if( list.size() == 1 ) {
				return new AlphabetContainer( alphabet[start[0]] );
			} else {
				return new AlphabetContainer( list.toArray( new Alphabet[0] ), ind );
			}
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.InstantiableFromParameterSet#getCurrentParameterSet()
	 */
	public AlphabetContainerParameterSet getCurrentParameterSet() throws Exception {
		if( parameters != null ) {
			return parameters.clone();
		} else {
			if( isSimple() ) {
				return new AlphabetContainerParameterSet( alphabet[0] );
			} else if( index.length == alphabet.length ) {
				return new AlphabetContainerParameterSet( alphabet );
			} else {
				return new AlphabetContainerParameterSet( alphabet, index );
			}
		}
	}

	/**
	 * Returns the delimiter that should be used (for writing e.g. a sequence).
	 * 
	 * @return the delimiter
	 */
	public String getDelim() {
		return delim;
	}

	/**
	 * Returns the maximal {@link Alphabet} length of this
	 * {@link AlphabetContainer}.
	 * 
	 * @return the maximal {@link Alphabet} length
	 */
	public double getMaximalAlphabetLength() {
		return l;
	}

	/**
	 * Returns the minimal value of the underlying {@link Alphabet} of position
	 * <code>pos</code>.
	 * 
	 * @param pos
	 *            the given position
	 * 
	 * @return the minimal value of the {@link Alphabet} of the given position
	 * 
	 * @see Alphabet#getMin()
	 */
	public double getMin( int pos ) {
		if( alphabet.length == 1 ) {
			return alphabet[0].getMin();
		} else {
			return alphabet[index[pos]].getMin();
		}
	}

	/**
	 * Returns the minimal {@link Alphabet} length of this
	 * {@link AlphabetContainer}.
	 * 
	 * @return the minimal {@link Alphabet} length of this
	 *         {@link AlphabetContainer}
	 */
	public double getMinimalAlphabetLength() {
		double length = alphabet[0].length();
		for( int i = 1; i < alphabet.length; i++ ) {
			length = Math.min( length, alphabet[i].length() );
		}
		return length;
	}

	/**
	 * Returns the possible length for {@link Sequence}s using this
	 * {@link AlphabetContainer}. If 0 (zero) is returned, all lengths are
	 * possible.
	 * 
	 * @return the possible length using this {@link AlphabetContainer}
	 */
	public int getPossibleLength() {
		return alphabet.length == 1 ? 0 : index.length;
	}

	/**
	 * Returns a sub-container with the {@link Alphabet}s for the positions
	 * starting at <code>start</code> and with length <code>length</code>. The
	 * method can be used for subsequences, ... .
	 * 
	 * @param start
	 *            the index of the start position
	 * @param length
	 *            the length
	 * 
	 * @return the sub-container of {@link Alphabet}s
	 * 
	 * @see AlphabetContainer#getCompositeContainer(int[], int[])
	 */
	public AlphabetContainer getSubContainer( int start, int length ) {
		if( alphabet.length == 1 || ( start == 0 && length == getPossibleLength() ) ) {
			return this;
		} else {
			int[] ind = new int[length];
			int[] used = new int[alphabet.length];
			Arrays.fill( used, -1 );
			ArrayList<Alphabet> list = new ArrayList<Alphabet>();
			for( int i = 0; i < length; i++, start++ ) {
				if( used[index[start]] < 0 ) {
					used[index[start]] = list.size();
					list.add( alphabet[index[start]] );
				}
				ind[i] = used[index[start]];
			}
			if( list.size() == 1 ) {
				return new AlphabetContainer( list.get( 0 ) );
			} else {
				return new AlphabetContainer( list.toArray( new Alphabet[0] ), ind );
			}
		}
	}

	/**
	 * Returns a {@link String} representation of the encoded symbol
	 * <code>val</code> of the {@link Alphabet} of position <code>pos</code> of
	 * this {@link AlphabetContainer}.
	 * 
	 * @param pos
	 *            the position of the {@link Alphabet}
	 * @param val
	 *            the value of the encoded symbol
	 * 
	 * @return a {@link String} representation for the encoded symbol
	 *         <code>val</code> of the {@link Alphabet} of position
	 *         <code>pos</code>
	 */
	public String getSymbol( int pos, double val ) {
		if( isDiscreteAt( pos ) ) {
			if( isSimple() ) {
				return ( (DiscreteAlphabet)alphabet[0] ).getSymbolAt( (int)val );
			} else {
				return ( (DiscreteAlphabet)alphabet[index[pos]] ).getSymbolAt( (int)val );
			}
		} else {
			return "" + val;
		}
	}

	/**
	 * Indicates if all used {@link Alphabet}s ignore the case.
	 * 
	 * @return <code>true</code> if all used alphabets ignore the case,
	 *         <code>false</code> otherwise
	 */
	public final boolean ignoresCase() {
		int i = 0;
		while( i < alphabet.length && ( alphabet[i] instanceof ContinuousAlphabet || ( (DiscreteAlphabet)alphabet[i] ).ignoresCase() ) ) {
			i++;
		}
		return alphabet.length == i;
	}

	/**
	 * Indicates if all positions use discrete {@link Alphabet}s.
	 * 
	 * @return <code>true</code> if all positions use discrete {@link Alphabet}s,
	 *         otherwise <code>false</code>
	 */
	public final boolean isDiscrete() {
		return getType() == AlphabetContainerType.DISCRETE;
	}

	/**
	 * Indicates if position <code>pos</code> is a discrete random variable,
	 * i.e. if the {@link Alphabet} of position <code>pos</code> is discrete.
	 * 
	 * @param pos
	 *            the position
	 * 
	 * @return <code>true</code> if position <code>pos</code> is a discrete
	 *         random variable, <code>false</code> otherwise
	 */
	public boolean isDiscreteAt( int pos ) {
		if( alphabet.length == 1 ) {
			return alphabet[0] instanceof DiscreteAlphabet;
		} else {
			return alphabet[index[pos]] instanceof DiscreteAlphabet;
		}
	}

	/**
	 * Indicates if <code>continuous</code> is a symbol of the {@link Alphabet}
	 * used at position <code>pos</code> of the {@link AlphabetContainer}.
	 * 
	 * @param pos
	 *            the position
	 * @param continuous
	 *            the continuous value
	 * 
	 * @return <code>true</code> if <code>continuous</code> is a symbol of the
	 *         {@link Alphabet} used in position <code>pos</code>,
	 *         <code>false</code> otherwise
	 */
	public boolean isEncodedSymbol( int pos, double continuous ) {
		if( isDiscreteAt( pos ) ) {
			int discrete = toDiscrete( pos, continuous );
			return ( (DiscreteAlphabet)getAlphabetAt( pos ) ).isEncodedSymbol( discrete ) && ( discrete - continuous == 0d );
		} else {
			return ( (ContinuousAlphabet)getAlphabetAt( pos ) ).isEncodedSymbol( continuous );
		}
	}

	/**
	 * Indicates whether all random variables are defined over the same range,
	 * i.e. if the {@link AlphabetContainer} is simple and all positions use the
	 * same (fixed) {@link Alphabet}.
	 * 
	 * @return whether all random variables are defined over the same range,
	 *         i.e. if the {@link AlphabetContainer} is simple
	 */
	public final boolean isSimple() {
		return alphabet.length == 1;
	}

	/**
	 * This method helps to determine if the {@link AlphabetContainer} also
	 * computes the reverse complement of a {@link Sequence}.
	 * 
	 * @return <code>true</code> if the {@link AlphabetContainer} also computes
	 *         the reverse complement of a {@link Sequence}, <code>false</code>
	 *         otherwise
	 */
	public final boolean isReverseComplementable() {
		return isSimple() && alphabet[0] instanceof ComplementableDiscreteAlphabet;
	}

	/**
	 * Returns the discrete value for <code>val</code> of the {@link Alphabet}
	 * of position <code>pos</code> in the {@link AlphabetContainer}.
	 * 
	 * @param pos
	 *            the position
	 * @param val
	 *            the value
	 * 
	 * @return the discrete value for <code>val</code> of the {@link Alphabet}
	 *         of position <code>pos</code>
	 * 
	 * @see AlphabetContainer#isDiscreteAt(int)
	 */
	public int toDiscrete( int pos, double val ) {
		if( isDiscreteAt( pos ) ) {
			return (int)val;
		} else {
			// TODO NaNs
			if(Double.isNaN( val )){
				throw new RuntimeException( "NaNs cannot be discretized" );
			}
			// TODO make it better (!?!, better discritisation)
			return (int)( val - getAlphabetAt( pos ).getMin() );
		}
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		String erg = "possible length: " + getPossibleLength() + "\n";
		erg += "alphabet : ";
		if( getPossibleLength() == 0 ) {
			erg += alphabet[0].toString();
		} else {
			for( int i = 0; i < getPossibleLength(); i++ ) {
				erg += "\n\t" + i + "\t" + getAlphabetAt( i ).toString();
			}
		}
		return erg + "\n";
	}

	/* (non-Javadoc)
	 * @see de.jstacs.Storable#toXML()
	 */
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer( 300 + alphabet.length * 200 );
		XMLParser.appendObjectWithTags( xml, alphabet, "Alphabets" );
		if( alphabet.length > 1 ) {
			XMLParser.appendObjectWithTags( xml, index, "Assignment" );
		}
		XMLParser.addTags( xml, XML_TAG );
		return xml;
	}

	/**
	 * Returns the type of this {@link AlphabetContainer}.
	 * 
	 * @return the type
	 * 
	 * @see AlphabetContainerType
	 */
	public final AlphabetContainerType getType() {
		return AlphabetContainerType.determineType( alphabet );
	}
	
	/**
	 * This method returns the index of the {@link Alphabet} that is used for the given position.
	 * 
	 * @param pos the position
	 * 
	 * @return the index of the used {@link Alphabet}
	 */
	public int getAlphabetIndexForPosition( int pos ) {
		if( alphabet.length == 1 ) {
			return 0;
		} else {
			return index[pos];
		}
	}
	
	/**
	 * This method returns the number of {@link Alphabet}s used in the current {@link AlphabetContainer}.
	 * 
	 * @return the number of used {@link Alphabet}s
	 */
	public int getNumberOfAlphabets() {
		return alphabet.length;
	}
	
	/**
	 * This method returns an object that is used for assigning the positions of the {@link Sequence}s to specific {@link Alphabet}s. 
	 * 
	 * @return <code>null</code> or an int array
	 */
	public int[] getIndexForAlphabets() {
		return index == null ? null : index.clone();
	}
}

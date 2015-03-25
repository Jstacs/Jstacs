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

import java.util.Iterator;
import java.util.LinkedList;

import de.jstacs.DataType;
import de.jstacs.data.AlphabetContainer.AbstractAlphabetContainerParameterSet;
import de.jstacs.data.AlphabetContainer.AlphabetContainerType;
import de.jstacs.data.alphabets.Alphabet;
import de.jstacs.data.alphabets.ContinuousAlphabet;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.alphabets.ContinuousAlphabet.ContinuousAlphabetParameterSet;
import de.jstacs.data.alphabets.DNAAlphabet.DNAAlphabetParameterSet;
import de.jstacs.data.alphabets.ProteinAlphabet.ProteinAlphabetParameterSet;
import de.jstacs.data.alphabets.DiscreteAlphabet.DiscreteAlphabetParameterSet;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.ArrayParameterSet;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.InstanceParameterSet;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;

/**
 * Class for the {@link ParameterSet} of an {@link AlphabetContainer}.
 * 
 * @author Jan Grau
 */
public class AlphabetContainerParameterSet extends AbstractAlphabetContainerParameterSet<AlphabetContainer> {

	private AlphabetContainerType type;

	private boolean simple;

	/**
	 * Creates a new {@link AlphabetContainerParameterSet} of an
	 * {@link AlphabetContainer} with {@link AlphabetContainerType}
	 * <code>type</code>. If <code>simple</code> is <code>true</code>, only a
	 * single {@link Alphabet} is expected.
	 * 
	 * @param type
	 *            the type of the alphabet(s)
	 * @param simple
	 *            indicates if there shall be only a single {@link Alphabet}
	 * @throws Exception if the parameters could not be created for the given arguments
	 * 
	 * @see AlphabetContainerType
	 * @see InstanceParameterSet#InstanceParameterSet(Class)
	 */
	public AlphabetContainerParameterSet( AlphabetContainerType type, boolean simple ) throws Exception {
		super( AlphabetContainer.class );
		this.type = type;
		this.simple = simple;
		addParameters();
	}
	
	/**
	 * Creates a new {@link AlphabetContainerParameterSet} of a
	 * complex {@link AlphabetContainer} with <code>type</code> {@link AlphabetContainerType#BOTH}.
	 * 
	 * @throws Exception if the {@link ParameterSet} could not be created
	 * 
	 * @see #AlphabetContainerParameterSet(AlphabetContainer.AlphabetContainerType, boolean)
	 */
	public AlphabetContainerParameterSet() throws Exception {
		this( AlphabetContainerType.BOTH, false );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link AlphabetContainerParameterSet} out of its XML
	 * representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link AlphabetContainerParameterSet} could not be
	 *             reconstructed out of the XML representation (the
	 *             {@link StringBuffer} <code>representation</code> could not be
	 *             parsed)
	 * 
	 * @see InstanceParameterSet#InstanceParameterSet(StringBuffer)
	 * @see de.jstacs.Storable
	 */
	public AlphabetContainerParameterSet( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}

	/**
	 * Creates a new {@link AlphabetContainerParameterSet} of a simple
	 * {@link AlphabetContainer} from a single {@link Alphabet}.
	 * 
	 * @param alph
	 *            the {@link Alphabet}
	 * 
	 * @throws Exception
	 *             if an error occurred during the creation of the appropriate
	 *             parameters
	 * 
	 * @see InstanceParameterSet#InstanceParameterSet(Class)
	 */
	public AlphabetContainerParameterSet( Alphabet alph ) throws Exception {
		super( AlphabetContainer.class );
		if( alph instanceof DiscreteAlphabet ) {
			this.type = AlphabetContainerType.DISCRETE;
		} else {
			this.type = AlphabetContainerType.CONTINUOUS;
		}
		this.simple = true;
		addParameters();
		if( type == AlphabetContainerType.CONTINUOUS ) {
			this.parameters.get( 0 ).setValue( alph.getCurrentParameterSet() );
		} else {
			InstanceParameterSet<? extends Alphabet> ap = alph.getCurrentParameterSet();
			this.parameters.get( 0 ).setValue( ap );
			( (SelectionParameter)this.parameters.get( 0 ) ).getParametersInCollection()
					.getParameterAt( ( (SelectionParameter)this.parameters.get( 0 ) ).getSelected() )
					.setValue( ap );
		}
	}

	/**
	 * Creates a new {@link AlphabetContainerParameterSet} from an array of
	 * {@link Alphabet}s.
	 * 
	 * @param alphabets
	 *            the array of {@link Alphabet}s
	 * 
	 * @throws Exception
	 *             if an error occurred during the creation of the appropriate
	 *             parameters
	 * 
	 * @see InstanceParameterSet#InstanceParameterSet(Class)
	 */
	public AlphabetContainerParameterSet( Alphabet[] alphabets ) throws Exception {
		this( AlphabetContainer.class, alphabets );
	}
	
	/**
	 * /**
	 * Creates a new {@link AlphabetContainerParameterSet} from an array of
	 * {@link Alphabet}s for a given sub-class of {@link AlphabetContainer}.
	 * @param instanceClass the sub-class fo {@link AlphabetContainer}
	 * @param alphabets the alphabets
	 * @throws Exception if the parameters could not be created for the given arguments
	 */
	protected AlphabetContainerParameterSet( Class<? extends AlphabetContainer> instanceClass, Alphabet... alphabets ) throws Exception {
		super( instanceClass );
		this.simple = false;
		this.type = AlphabetContainerType.determineType( alphabets );
		addParameters();
		AlphabetArrayParameterSet pars = new AlphabetArrayParameterSet( alphabets, this.type );
		parameters.get( 0 ).setValue( pars );
		( (SelectionParameter)parameters.get( 0 ) ).getParametersInCollection()
				.getParameterAt( ( (SelectionParameter)parameters.get( 0 ) ).getSelected() )
				.setValue( pars );
	}

	/**
	 * Creates a new {@link AlphabetContainerParameterSet} from an array of
	 * {@link Alphabet}s and an array of <code>int</code>s defining the
	 * {@link Alphabet} index <code>i</code> in <code>alphabets</code>
	 * that is used for position <code>i</code>.
	 * 
	 * @param alphabets
	 *            the {@link Alphabet}s
	 * @param indexes
	 *            the indexes
	 * 
	 * @throws Exception
	 *             if an error occurred during the creation of the appropriate
	 *             parameters
	 * 
	 * @see InstanceParameterSet#InstanceParameterSet(Class)
	 */
	public AlphabetContainerParameterSet( Alphabet[] alphabets, int[] indexes ) throws Exception {
		super( AlphabetContainer.class );
		this.simple = false;
		this.type = AlphabetContainerType.determineType( alphabets );
		addParameters();
		SectionDefinedAlphabetParameterSet pars = new SectionDefinedAlphabetParameterSet( alphabets, indexes );
		parameters.get( 0 ).setValue( pars );
		( (SelectionParameter)parameters.get( 0 ) ).getParametersInCollection()
				.getParameterAt( ( (SelectionParameter)parameters.get( 0 ) ).getSelected() )
				.setValue( pars );
	}

	/**
	 * Indicates if all positions use {@link DiscreteAlphabetParameterSet}, i.e.
	 * if all {@link Alphabet}s of the corresponding {@link AlphabetContainer}
	 * are discrete.
	 * 
	 * @return <code>true</code> if all positions are discrete,
	 *         <code>false</code> otherwise
	 */
	public boolean isDiscrete() {
		return type == AlphabetContainerType.DISCRETE;
	}

	public boolean isSimple() {
		return simple;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.parameters.InstanceParameterSet#toXML()
	 */
	@Override
	public StringBuffer toXML() {
		StringBuffer buf = super.toXML();
		XMLParser.addTags( buf, "superParameters" );
		XMLParser.appendObjectWithTags( buf, type, "type" );
		XMLParser.appendObjectWithTags( buf, simple, "simple" );
		XMLParser.addTags( buf, "alphabetContainerParameterSet" );
		return buf;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.parameters.InstanceParameterSet#fromXML(java.lang.StringBuffer)
	 */
	@Override
	protected void fromXML( StringBuffer representation ) throws NonParsableException {
		representation = XMLParser.extractForTag( representation, "alphabetContainerParameterSet" );
		super.fromXML( XMLParser.extractForTag( representation, "superParameters" ) );
		type = XMLParser.extractObjectForTags( representation, "type", AlphabetContainerType.class );
		simple = XMLParser.extractObjectForTags( representation, "simple", boolean.class );
	}

	public int getPossibleLength() {
		Object o = parameters.get( 0 ).getValue();
		if( o instanceof ProteinAlphabetParameterSet || o instanceof DNAAlphabetParameterSet || o instanceof DiscreteAlphabetParameterSet
			|| o instanceof ContinuousAlphabetParameterSet ) {
			return 0;
		} else {
			ParameterSet set = (ParameterSet)o;
			return (Integer)set.getParameterAt( 0 ).getValue();
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.parameters.ParameterSet#loadParameters()
	 */
	private void addParameters() throws Exception {
		LinkedList<InstanceParameterSet> list = type.getInstanceParameterSets();
		
		ParameterSet[] ar = new ParameterSet[list.size()+(simple?0:2)];
		list.toArray( ar );
		
		if( !simple ) {
			AlphabetArrayParameterSet arrayParameters = new AlphabetArrayParameterSet( type );
			SectionDefinedAlphabetParameterSet sectionParameters = new SectionDefinedAlphabetParameterSet( type );
			
			ar[list.size()] = arrayParameters;
			ar[list.size()+1] = sectionParameters;
		}

		parameters.add( new SelectionParameter( "Alphabet",
					simple ? (
							type==AlphabetContainerType.DISCRETE ? "Select a discrete alphabet" :
								type==AlphabetContainerType.CONTINUOUS ? "Set the continuous alphabet" :
									"Select an alphabet"
					) : (
							type==AlphabetContainerType.DISCRETE ? "Select the discrete alphabets" :
								type==AlphabetContainerType.CONTINUOUS ? "Set the continuous alphabets" :
									"Select the alphabets"
					),
					true, ar ) );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.parameters.InstanceParameterSet#getInstanceName()
	 */
	@Override
	public String getInstanceName() {
		return "Alphabet";
	}

	/* (non-Javadoc)
	 * @see de.jstacs.parameters.InstanceParameterSet#getInstanceComment()
	 */
	@Override
	public String getInstanceComment() {
		return "Set an alphabet for your model or data.";
	}

	/**
	 * Class for the parameter set of an array of {@link Alphabet}s where each
	 * {@link Alphabet} may be used for one or more sections of positions.
	 * 
	 * @author Jan Grau
	 * 
	 */
	public static class SectionDefinedAlphabetParameterSet extends ExpandableParameterSet {

		private AlphabetContainerType type;

		/**
		 * Creates a new {@link SectionDefinedAlphabetParameterSet} for a set of
		 * discrete or continuous {@link Alphabet}s.
		 * 
		 * @param type
		 *            the type of the {@link Alphabet}(s)
		 * 
		 * @throws Exception
		 *             if the {@link SectionDefinedAlphabetParameterSet} could
		 *             not be created
		 * 
		 * @see ExpandableParameterSet#ExpandableParameterSet(ParameterSet,
		 *      String, String)
		 */
		public SectionDefinedAlphabetParameterSet( AlphabetContainerType type ) throws Exception {
			super( new SimpleParameterSet( new Parameter[]{
					(type == AlphabetContainerType.CONTINUOUS
						? //continuous
							new ParameterSetContainer( "Alphabet",
										"Set the parameters of the alphabet.",
										new ContinuousAlphabet.ContinuousAlphabetParameterSet() )
						: //discrete || both
							new SelectionParameter( "Type of alphabet", "Select the type of the alphabet",
										true, type.getInstanceParameterSets().toArray(new InstanceParameterSet[0]) )
							),
					new SimpleParameter( DataType.STRING,
										"Section",
										"The section for that this alphabet is applied. Use &quot;,&quot; to separate different positions and &quot;-&quot; to indicate ranges, e.g. 1-3,5.",
										true ) } ),
				"Alphabet",
				"Set the alphabet", 0 );
			SimpleParameter length = new SimpleParameter( DataType.INT, "Length", "The length of the array.", true );
			length.setRangeable( false );
			this.parameters.add( 0, length );
			this.type = type;
		}
		
		/**
		 * Creates a new {@link SectionDefinedAlphabetParameterSet} for a set of
		 * discrete and continuous {@link Alphabet}s.
		 * 
		 * @throws Exception
		 *             if the {@link SectionDefinedAlphabetParameterSet} could
		 *             not be created
		 */
		public SectionDefinedAlphabetParameterSet() throws Exception {
			this( AlphabetContainerType.BOTH );
		}

		/**
		 * Creates a new {@link SectionDefinedAlphabetParameterSet} from an
		 * array of {@link Alphabet}s and an array of indexes that define the
		 * index of the {@link Alphabet} in <code>alphabets</code> belonging to
		 * that position in <code>indexes</code>.
		 * 
		 * @param alphabets
		 *            the array of {@link Alphabet}s
		 * @param indexes
		 *            the array of indexes of the {@link Alphabet}s belonging to
		 *            the positions
		 * 
		 * @throws Exception
		 *             if the {@link SectionDefinedAlphabetParameterSet} could
		 *             not be created
		 * 
		 * @see de.jstacs.data.AlphabetContainerParameterSet.SectionDefinedAlphabetParameterSet#AlphabetContainerParameterSet.SectionDefinedAlphabetParameterSet(AlphabetContainer.AlphabetContainerType) 
		 */
		public SectionDefinedAlphabetParameterSet( Alphabet[] alphabets, int[] indexes ) throws Exception {
			this( AlphabetContainerType.determineType( alphabets ) );
			parameters.get( 0 ).setValue( indexes.length );
			for( int i = 0; i < alphabets.length; i++ ) {
				addParameterToSet();
			}

			for( int i = 0; i < alphabets.length; i++ ) {

				ParameterSet temp = (ParameterSet)parameters.get( i + 1 ).getValue();
				if( type == AlphabetContainerType.CONTINUOUS ) {
					temp.getParameterAt( 0 ).setValue( alphabets[i].getCurrentParameterSet() );
				} else {
					InstanceParameterSet ap = alphabets[i].getCurrentParameterSet();
					temp.getParameterAt( 0 ).setValue( ap );
					( (SelectionParameter)temp.getParameterAt( 0 ) ).getParametersInCollection()
							.getParameterAt( ( (SelectionParameter)temp.getParameterAt( 0 ) ).getSelected() )
							.setValue( ap );
				}
				StringBuffer sections = new StringBuffer();
				for( int j = 0; j < indexes.length; j++ ) {
					if( indexes[j] == i ) {
						sections.append( ( j + 1 ) + ", " );
					}
				}
				sections.delete( sections.length() - 2, sections.length() );
				temp.getParameterAt( 1 ).setValue( sections.toString() );
			}
		}

		/**
		 * The standard constructor for the interface {@link de.jstacs.Storable}
		 * . Creates a new {@link SectionDefinedAlphabetParameterSet} out of its
		 * XML representation.
		 * 
		 * @param representation
		 *            the XML representation as {@link StringBuffer}
		 * 
		 * @throws NonParsableException
		 *             if the
		 *             {@link de.jstacs.data.alphabets.Alphabet.AlphabetParameterSet}
		 *             could not be reconstructed out of the XML representation
		 *             (the {@link StringBuffer} <code>representation</code>
		 *             could not be parsed)
		 * 
		 * @see ExpandableParameterSet#ExpandableParameterSet(StringBuffer)
		 * @see de.jstacs.Storable
		 */
		public SectionDefinedAlphabetParameterSet( StringBuffer representation ) throws NonParsableException {
			super( representation );
		}

		/* (non-Javadoc)
		 * @see de.jstacs.parameters.ExpandableParameterSet#toXML()
		 */
		@Override
		public StringBuffer toXML() {
			StringBuffer buf = super.toXML();
			XMLParser.addTags( buf, "superParameters" );
			XMLParser.appendObjectWithTags( buf, type, "type" );
			XMLParser.addTags( buf, "sectionDefinedAlphabetParameterSet" );
			return buf;
		}

		/* (non-Javadoc)
		 * @see de.jstacs.parameters.ExpandableParameterSet#fromXML(java.lang.StringBuffer)
		 */
		@Override
		public void fromXML( StringBuffer representation ) throws NonParsableException {
			representation = XMLParser.extractForTag( representation, "sectionDefinedAlphabetParameterSet" );
			super.fromXML( XMLParser.extractForTag( representation, "superParameters" ) );
			type = XMLParser.extractObjectForTags( representation, "type", AlphabetContainerType.class );
		}

		/* (non-Javadoc)
		 * @see de.jstacs.parameters.ParameterSet#hasDefaultOrIsSet()
		 */
		@Override
		public boolean hasDefaultOrIsSet() {
			if( !super.hasDefaultOrIsSet() ) {
				return false;
			} else {
				boolean[] set = new boolean[(Integer)parameters.get( 0 ).getValue()];
				LinkedList<Integer> tempList;
				int no = 0;
				for(int i=1;i<parameters.size();i++){
					Parameter element = parameters.get( i );
					no++;
					ParameterSet temp = (ParameterSet)element.getValue();
					String section = (String)temp.getParameterAt( 1 ).getValue();
					try {
						tempList = parseSections( section );
						Iterator<Integer> posIt = tempList.iterator();
						while( posIt.hasNext() ) {
							int curr = posIt.next();
							if( curr >= set.length ) {
								errorMessage = "Position " + ( curr + 1 ) + " out of range defined by length.";
								return false;
							} else if( set[curr] ) {
								errorMessage = "Alphabet for position " + ( curr + 1 ) + " defined at least twice.";
								return false;
							} else {
								set[curr] = true;
							}

						}

					} catch ( Exception e ) {
						errorMessage = "Malformed section definition no. " + no;
						e.printStackTrace();
						return false;
					}
				}
				for( int i = 0; i < set.length; i++ ) {
					if( !set[i] ) {
						errorMessage = "No alphabet defined for position " + ( i + 1 ) + ".";
						return false;
					}
				}
				errorMessage = null;
				return true;
			}
		}

		/**
		 * Parses the sections as defined in <code>sections</code> to a
		 * {@link LinkedList} of <code>Integer</code>s. Sections may be defined
		 * as ranges (e.g. &quot;3-6&quot;), as lists (e.g. &quot;3,4,5,6&quot;)
		 * or as combinations of both (e.g. &quot;3-5,6&quot;).
		 * 
		 * @param sections
		 *            the sections
		 * 
		 * @return the {@link LinkedList} of positions
		 * 
		 * @throws Exception
		 *             if <code>sections</code> could not be parsed, contains a
		 *             position more than once or contains positions that are
		 *             out of range
		 */
		public static LinkedList<Integer> parseSections( String sections ) throws Exception {
			String[] secs = sections.split( "(\\s*,\\s*)" );
			LinkedList<Integer> list = new LinkedList<Integer>();
			String[] temp;
			for( int i = 0; i < secs.length; i++ ) {
				// if(secs[i].indexOf('-') > -1){
				temp = secs[i].split( "(\\s*-\\s*)" );
				// }else{
				// temp = new String[]{secs[i]};
				// }
				if( temp.length == 2 ) {
					int start = Integer.parseInt( temp[0] ) - 1;
					int end = Integer.parseInt( temp[1] ) - 1;
					start = start <= end ? start : end;
					end = end >= start ? end : start;
					if( start < 0 ) {
						throw new Exception( "Malformed sections, no negative values allowed." );
					}
					for( int j = start; j <= end; j++ ) {
						list.add( j );
					}
				} else if( temp.length == 1 ) {
					int pos = Integer.parseInt( temp[0] ) - 1;
					if( pos < 0 ) {
						throw new Exception( "Malformed sections, no negative values allowed." );
					}
					list.add( pos );
				} else {
					throw new Exception( "Malformed sections." );
				}
			}
			return list;
		}

		/* (non-Javadoc)
		 * @see de.jstacs.parameters.ExpandableParameterSet#clone()
		 */
		@Override
		public SectionDefinedAlphabetParameterSet clone() throws CloneNotSupportedException {
			try {
				SectionDefinedAlphabetParameterSet clone = (SectionDefinedAlphabetParameterSet)super.clone();
				
				return clone;
			} catch ( Exception e ) {
				e.printStackTrace();
				throw new CloneNotSupportedException( e.getCause().getMessage() );
			}

		}

		/**
		 * Returns a descriptive name for this
		 * {@link SectionDefinedAlphabetParameterSet}.
		 * 
		 * @return the descriptive name
		 */
		public String getInstanceName() {
			return "Alphabets defined by section";
		}

		/**
		 * Returns a descriptive comment on this
		 * {@link SectionDefinedAlphabetParameterSet}.
		 * 
		 * @return the descriptive comment
		 */
		public String getInstanceComment() {
			return "Set the alphabets for all positions.";
		}

	}

	/**
	 * Class for the parameters of an array of {@link Alphabet}s of defined
	 * length.
	 * 
	 * @author Jan Grau
	 */
	public static class AlphabetArrayParameterSet extends ArrayParameterSet {

		private AlphabetContainerType type;

		/**
		 * Creates a new {@link AlphabetArrayParameterSet} from the information
		 * about the <code>type</code> of the {@link Alphabet}s, e.g. if the
		 * array shall contain only the parameters for discrete {@link Alphabet}
		 * s.
		 * 
		 * @param type
		 *            the type of the {@link Alphabet}(s)
		 * 
		 * @throws Exception
		 *             if the {@link AlphabetArrayParameterSet} could not be
		 *             created
		 * 
		 * @see ArrayParameterSet#ArrayParameterSet(ParameterSet, String,
		 *      String)
		 */
		public AlphabetArrayParameterSet( AlphabetContainerType type ) throws Exception {
			super( (
			
			type == AlphabetContainerType.CONTINUOUS
					? //continuous
							new ContinuousAlphabet.ContinuousAlphabetParameterSet()
					: //discrete || both
						new SimpleParameterSet( new Parameter[]{ new SelectionParameter( "Type of alphabet", "Select the type of the alphabet",
																						true, type.getInstanceParameterSets().toArray(new InstanceParameterSet[0]) )
												} )

			),
					"Alphabet",
					"Set the alphabet" );
		}
		
		/**
		 * Creates a new {@link AlphabetArrayParameterSet} with <code>type</code> {@link de.jstacs.data.AlphabetContainer.AlphabetContainerType#BOTH}.
		 * 
		 * @throws Exception
		 *             if the {@link AlphabetArrayParameterSet} could not be
		 *             created
		 */
		public AlphabetArrayParameterSet() throws Exception {
			this( AlphabetContainerType.BOTH );
		}

		/**
		 * The standard constructor for the interface {@link de.jstacs.Storable}
		 * . Creates a new {@link AlphabetArrayParameterSet} out of its XML
		 * representation.
		 * 
		 * @param representation
		 *            the XML representation as {@link StringBuffer}
		 * 
		 * @throws NonParsableException
		 *             if the {@link AlphabetArrayParameterSet} could not be
		 *             reconstructed out of the XML representation (the
		 *             {@link StringBuffer} <code>representation</code> could
		 *             not be parsed)
		 * 
		 * @see ArrayParameterSet#ArrayParameterSet(StringBuffer)
		 * @see de.jstacs.Storable
		 */
		public AlphabetArrayParameterSet( StringBuffer representation ) throws NonParsableException {
			super( representation );
		}

		/**
		 * Creates a new {@link AlphabetArrayParameterSet} from an array of
		 * {@link Alphabet}s and the information about the <code>type</code> of
		 * the {@link Alphabet}s.
		 * 
		 * @param alphabets
		 *            the {@link Alphabet}s
		 * @param type
		 *            the type of the {@link AlphabetContainer}
		 * 
		 * @throws Exception
		 *             if the {@link AlphabetArrayParameterSet} could not be
		 *             created
		 */
		public AlphabetArrayParameterSet( Alphabet[] alphabets, AlphabetContainerType type ) throws Exception {
			this( type );
			this.parameters.get( 0 ).setValue( alphabets.length );
			for( int i = 0; i < alphabets.length; i++ ) {
				if( type == AlphabetContainerType.CONTINUOUS ) {
					parameters.get( i + 1 ).setValue( alphabets[i].getCurrentParameterSet() );
				} else {
					InstanceParameterSet ap = alphabets[i].getCurrentParameterSet();
					( (ParameterSet)parameters.get( i + 1 ).getValue() ).getParameterAt( 0 ).setValue( ap );
					( (SelectionParameter)( (ParameterSet)parameters.get( i + 1 ).getValue() ).getParameterAt( 0 ) ).getParametersInCollection()
							.getParameterAt( ( (SelectionParameter)( (ParameterSet)parameters.get( i + 1 ).getValue() ).getParameterAt( 0 ) ).getSelected() )
							.setValue( ap );
				}
			}
		}

		/* (non-Javadoc)
		 * @see de.jstacs.parameters.ExpandableParameterSet#clone()
		 */
		@Override
		public AlphabetArrayParameterSet clone() throws CloneNotSupportedException {
			try {
				AlphabetArrayParameterSet clone = (AlphabetArrayParameterSet)super.clone();
				return clone;
			} catch ( Exception e ) {
				throw new CloneNotSupportedException( e.getCause().getMessage() );
			}
		}

		/* (non-Javadoc)
		 * @see de.jstacs.parameters.ArrayParameterSet#toXML()
		 */
		@Override
		public StringBuffer toXML() {
			StringBuffer buf = super.toXML();
			XMLParser.addTags( buf, "superParameters" );
			XMLParser.appendObjectWithTags( buf, type, "type" );
			XMLParser.addTags( buf, "alphabetArrayParameterSet" );
			return buf;
		}

		/* (non-Javadoc)
		 * @see de.jstacs.parameters.ArrayParameterSet#fromXML(java.lang.StringBuffer)
		 */
		@Override
		public void fromXML( StringBuffer representation ) throws NonParsableException {
			representation = XMLParser.extractForTag( representation, "alphabetArrayParameterSet" );
			super.fromXML( XMLParser.extractForTag( representation, "superParameters" ) );
			type = XMLParser.extractObjectForTags( representation, "type", AlphabetContainerType.class );
		}

		/**
		 * Returns a descriptive name for this {@link AlphabetArrayParameterSet}
		 * .
		 * 
		 * @return the descriptive name
		 */
		public String getInstanceName() {
			return "Alphabet-array";
		}

		/**
		 * Returns a descriptive comment on this
		 * {@link AlphabetArrayParameterSet}.
		 * 
		 * @return the descriptive comment
		 */
		public String getInstanceComment() {
			return "An array of alphabets where each position can have its own alphabet.";
		}
	}
}

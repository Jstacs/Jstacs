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

package de.jstacs.data.bioJava;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.SimpleAnnotation;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.io.NameTokenization;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.PointLocation;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.bio.symbol.SimpleAlphabet;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojavax.SimpleNamespace;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.io.RichSequenceBuilder;
import org.biojavax.bio.seq.io.RichSequenceBuilderFactory;
import org.biojavax.ontology.SimpleComparableTerm;

import de.jstacs.WrongAlphabetException;
import de.jstacs.data.Alphabet;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.Sequence;
import de.jstacs.data.alphabets.DNAAlphabet;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.annotation.LocatedSequenceAnnotation;
import de.jstacs.data.sequences.annotation.LocatedSequenceAnnotationWithLength;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.data.sequences.annotation.StrandedLocatedSequenceAnnotationWithLength;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.Result;

/**
 * This class provides static methods to convert BioJava datatypes (
 * {@link SequenceIterator}, {@link org.biojava.bio.seq.Sequence}) to
 * {@link DataSet}s and vice versa.
 * 
 * <br>
 * <br>
 * 
 * Here are two small examples how to create a
 * {@link org.biojavax.bio.seq.RichSequenceIterator} that can be used to create
 * a {@link DataSet}.
 * 
 * <br>
 * <br>
 * 
 * <code>
 * GenbankRichSequenceDB db = new GenbankRichSequenceDB();<br>
 * org.biojava.bio.seq.Sequence seq = db.getSequence( id );<br>
 * RichSequenceIterator iter = new RichSequence.IOTools.SingleRichSeqIterator(seq);<br>
 * </code> <br>
 * or <br>
 * <br>
 * <code>
 * RichSequenceIterator iter = RichSequence.IOTools.readFile( new File( fName ), null );
 * </code>
 * 
 * @author Jens Keilwagen, Jan Grau
 */
public class BioJavaAdapter {

	private static final Result[] EMPTY_RESULT_ARRAY = new Result[0];

	private static final String ANNOTATION_ID = "BJRSA";

	/**
	 * This method creates a new {@link DataSet} from a {@link SequenceIterator}.
	 * In cases where BioJava does return {@link org.biojava.bio.seq.Sequence}s
	 * instead of a {@link SequenceIterator}, you can use a
	 * {@link SimpleSequenceIterator} to wrap them up.
	 * 
	 * @param it
	 *            the sequence iterator
	 * @param filter
	 *            <code>null</code> or an arbitrary feature filter that
	 *            determines which features will be adopted
	 * 
	 * @return the {@link org.biojava.bio.seq.Sequence}s in <code>it</code>
	 *         converted to a {@link DataSet}
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public static DataSet sequenceIteratorToSample( SequenceIterator it, FeatureFilter filter ) throws Exception {
		LinkedList<Sequence> parsed = new LinkedList<Sequence>();
		org.biojava.bio.seq.Sequence current = it.nextSequence();
		org.biojava.bio.symbol.Alphabet abcCurrent, bioAbc = current.getAlphabet();
		List l = bioAbc.getAlphabets();
		Alphabet[] abc = new Alphabet[l.size()];
		for( int i = 0; i < abc.length; i++ ) {
			abcCurrent = (org.biojava.bio.symbol.Alphabet)l.get( i );
			if( abcCurrent.getAlphabets().size() == 1 ) {
				if( abcCurrent.getName().equals( "DNA" ) ) {
					abc[i] = new DNAAlphabet();
				} else if( abcCurrent instanceof FiniteAlphabet ) {
					FiniteAlphabet finAbc = (FiniteAlphabet)bioAbc;
					SymbolTokenization t = finAbc.getTokenization( "default" );
					Iterator symIt = finAbc.iterator();
					LinkedList<String> sym = new LinkedList<String>();
					Symbol s;
					while( symIt.hasNext() ) {
						s = (Symbol)symIt.next();
						//System.out.println( t.tokenizeSymbol(s) + "\t" + s.getAnnotation().toString() + "\t" + s.getName() + "\t" + s.toString() );
						sym.add( t.tokenizeSymbol( s ) );
					}
					abc[i] = new DiscreteAlphabet( true, sym.toArray( new String[0] ) );
				} else {
					//TODO
					throw new Exception( "not finite" );
				}
			} else {
				//TODO
				throw new Exception( "size > 1" );
			}
		}
		AlphabetContainer con = new AlphabetContainer( abc );

		SequenceAnnotation[] empty = new SequenceAnnotation[0];
		LinkedList<SequenceAnnotation> annotations = new LinkedList<SequenceAnnotation>();
		do {
			annotations.clear();
			if( current instanceof RichSequence ) {
				RichSequence rseq = (RichSequence)current;

				annotations.add( getRichSequenceAnnotation( rseq ) );
			}

			FeatureHolder features = null;
			if( filter == null ) {
				features = current;
			} else {
				features = current.filter( filter );
			}

			if( features.countFeatures() > 0 ) {
				Iterator<Feature> itFeat = features.features();
				while( itFeat.hasNext() ) {
					annotations.add( featuresToAnnotation( itFeat.next() ) );
				}
			}
			SequenceAnnotation[] annota = null;
			if( annotations.size() > 0 ) {
				annota = annotations.toArray( empty );
			}

			//TODO ambiguous symbols?
			parsed.add( Sequence.create( con, annota, current.seqString(), con.getDelim() ) );
			if( it.hasNext() ) {
				current = it.nextSequence();
			} else {
				break;
			}
		} while( true );
		return new DataSet( "", parsed.toArray( new Sequence[0] ) );
	}

	/**
	 * Returns the {@link SequenceAnnotation} corresponding to the properties
	 * that are associated with the {@link RichSequence} <code>rseq</code> via
	 * the methods {@link RichSequence#getDescription()},
	 * {@link RichSequence#getAccession()}, {@link RichSequence#getName()},
	 * {@link RichSequence#getIdentifier()} and
	 * {@link RichSequence#getVersion()}.
	 * 
	 * @param rseq
	 *            the sequence
	 * 
	 * @return the annotation of the sequence
	 * 
	 * @see SequenceAnnotation
	 */
	private static SequenceAnnotation getRichSequenceAnnotation( RichSequence rseq ) {
		LinkedList<Result> results = new LinkedList<Result>();
		String temp = rseq.getDescription();
		if( temp != null ) {
			results.add( new CategoricalResult( "Description", "The description of the sequence", temp ) );
		}
		temp = rseq.getAccession();
		if( temp != null ) {
			results.add( new CategoricalResult( "Accession", "The accession of the sequence", temp ) );
		}
		temp = rseq.getName();
		if( temp != null ) {
			results.add( new CategoricalResult( "Name", "The name of the sequence", temp ) );
		}
		temp = rseq.getIdentifier();
		if( temp != null ) {
			results.add( new CategoricalResult( "ID", "The identifier of the sequence", temp ) );
		}
		int vers = rseq.getVersion();
		results.add( new NumericalResult( "Version", "The version of the sequence", vers ) );

		return new SequenceAnnotation( "BioJava RichSequence Annotation", ANNOTATION_ID, results );
	}

	/**
	 * Converts <code>feature</code> and all of its sub-features to
	 * corresponding {@link SequenceAnnotation}s.
	 * 
	 * @param feature
	 *            the feature
	 * 
	 * @return the annotation object
	 * 
	 * @see SequenceAnnotation
	 */
	@SuppressWarnings( "unchecked" )
	private static SequenceAnnotation featuresToAnnotation( Feature feature ) {
		LinkedList<SequenceAnnotation> subAnnot = new LinkedList<SequenceAnnotation>();
		if( feature.countFeatures() > 0 ) {
			Iterator<Feature> it = feature.features();
			while( it.hasNext() ) {
				subAnnot.add( featuresToAnnotation( it.next() ) );
			}
		}

		String id = null;

		LinkedList<Result> results = new LinkedList<Result>();
		Annotation ann = feature.getAnnotation();
		Iterator<Map.Entry> it = ann.asMap().entrySet().iterator();
		while( it.hasNext() ) {
			Map.Entry en = it.next();
			if( en.getKey() instanceof SimpleComparableTerm ) {
				if( ( (SimpleComparableTerm)en.getKey() ).getName().equalsIgnoreCase( "identifier" ) ) {
					id = en.getValue().toString();
				} else {
					results.add( new CategoricalResult( ( (SimpleComparableTerm)en.getKey() ).getName(),
							( (SimpleComparableTerm)en.getKey() ).getDescription(),
							en.getValue().toString() ) );
				}
			} else {
				if( en.getKey().toString().equalsIgnoreCase( "identifier" ) ) {
					id = en.getValue().toString();
				} else {
					results.add( new CategoricalResult( en.getKey().toString(), "", en.getValue().toString() ) );
				}
			}
		}

		SequenceAnnotation[] subs = null;
		if( subAnnot.size() > 0 ) {
			subs = subAnnot.toArray( new SequenceAnnotation[0] );
		}

		Location loc = feature.getLocation();
		String type = feature.getType();
		if( id == null ) {
			id = feature.getType() + "_[" + feature.getLocation().toString() + "]";
		}
		org.biojava.bio.seq.Sequence fseq = feature.getSequence();
		if( fseq instanceof RichSequence ) {
			results.add( new CategoricalResult( "Accession", "Accession", ( (RichSequence)fseq ).getAccession() ) );
		}
		if( loc == Location.empty ) {
			return new SequenceAnnotation( type, id, subs, results.toArray( EMPTY_RESULT_ARRAY ) );
		} else if( loc.getMin() == loc.getMax() ) {
			return new LocatedSequenceAnnotation( loc.getMin() - 1, type, id, subs, results.toArray( EMPTY_RESULT_ARRAY ) );
		} else if( loc.isContiguous() ) {
			return new LocatedSequenceAnnotationWithLength( loc.getMin() - 1,
					loc.getMax() - loc.getMin() + 1,
					type,
					id,
					subs,
					results.toArray( EMPTY_RESULT_ARRAY ) );
		} else {
			Iterator<Location> subLocs = loc.blockIterator();
			LinkedList<SequenceAnnotation> subLocAnnots = new LinkedList<SequenceAnnotation>();
			int i = 0;
			while( subLocs.hasNext() ) {
				Location subLoc = subLocs.next();
				if( subLoc.getMin() == subLoc.getMax() ) {
					subLocAnnots.add( new LocatedSequenceAnnotation( subLoc.getMin() - 1, feature.getType(), "Sub location " + ( i + 1 ) ) );
				} else {
					subLocAnnots.add( new LocatedSequenceAnnotationWithLength( subLoc.getMin() - 1,
							subLoc.getMax() - subLoc.getMin() + 1,
							feature.getType(),
							"Sub location " + ( i + 1 ) ) );
				}
				i++;
			}
			subLocAnnots.addAll( 0, subAnnot );
			return new LocatedSequenceAnnotationWithLength( loc.getMin() - 1,
					loc.getMax() - loc.getMin() + 1,
					type,
					id,
					subLocAnnots.toArray( new SequenceAnnotation[0] ) );
		}

	}

	/**
	 * Creates a {@link SequenceIterator} from the {@link DataSet}
	 * <code>sample</code> preserving as much annotation as possible. This
	 * method works only for discrete alphabets.
	 * 
	 * @param sample
	 *            the {@link DataSet}
	 * @param flat
	 *            indicates if the features should be flattened, this may be
	 *            necessary to preserve all features, because some data formats
	 *            do not support hierarchical features
	 * 
	 * @return the corresponding {@link SequenceIterator}
	 * 
	 * @throws WrongAlphabetException
	 *             if the alphabet of the {@link DataSet} is not discrete
	 * @throws BioException
	 *             if forwarded from BioJava
	 */
	public static SequenceIterator sampleToSequenceIterator( DataSet sample, boolean flat ) throws WrongAlphabetException, BioException {
		LinkedList<org.biojava.bio.seq.Sequence> bjSeqs = new LinkedList<org.biojava.bio.seq.Sequence>();

		AlphabetContainer cont = sample.getAlphabetContainer();

		if( !( cont.isSimple() && cont.isDiscrete() ) ) {
			throw new WrongAlphabetException( "Only simple and discrete alphabets can be used" );
		}

		FiniteAlphabet bjAlph = null;

		if( cont.isSimple() ) {
			DiscreteAlphabet alph = (DiscreteAlphabet)cont.getAlphabetAt( 0 );
			if( alph instanceof DNAAlphabet ) {
				bjAlph = DNATools.getDNA();
			} else {
				Set<Symbol> symbols = new HashSet<Symbol>();
				for( int i = 0; i < (int)alph.length(); i++ ) {
					symbols.add( AlphabetManager.createSymbol( alph.getSymbolAt( i ) ) );
				}
				bjAlph = new SimpleAlphabet( symbols, "discrete" );
				AlphabetManager.registerAlphabet( bjAlph.getName(), bjAlph );
				( (SimpleAlphabet)bjAlph ).putTokenization( "default", new WorkingWordTokenization( bjAlph,
						!cont.ignoresCase(),
						cont.getDelim() ) );
			}
		} else {
			//TODO this part of the method will never be evaluated because of the "if" at the beginning of the method

			//We don't know how to create the homolog of an AlphabetContainer with more than 1 alphabet in BioJava
			//This is code might be helpfull to solve the problem.
			LinkedList<org.biojava.bio.symbol.Alphabet> alphabets = new LinkedList<org.biojava.bio.symbol.Alphabet>();
			for( int i = 0; i < cont.getPossibleLength(); i++ ) {
				DiscreteAlphabet alph = (DiscreteAlphabet)cont.getAlphabetAt( i );
				org.biojava.bio.symbol.Alphabet tempAlph = null;
				if( alph instanceof DNAAlphabet ) {
					tempAlph = DNATools.getDNA();
				} else {
					Set<Symbol> symbols = new HashSet<Symbol>();
					for( int j = 0; j < (int)alph.length(); j++ ) {
						symbols.add( AlphabetManager.createSymbol( alph.getSymbolAt( j ) ) );
					}
					tempAlph = new SimpleAlphabet( symbols );
				}
				alphabets.add( tempAlph );
			}
			//bjAlph = AlphabetManager.getCrossProductAlphabet( alphabets );
		}

		int numSeqs = sample.getNumberOfElements();
		for( int i = 0; i < numSeqs; i++ ) {
			Sequence seq = sample.getElementAt( i );

			RichSequence bjSeq = getSequence( bjAlph, cont, seq, "Sequence_" + ( i + 1 ) );
			SequenceAnnotation[] ann = seq.getAnnotation();
			if( ann != null ) {
				for( int a = 0; a < ann.length; a++ ) {
					if( !ann[a].getIdentifier().equals( ANNOTATION_ID ) ) {
						annotationToFeatures( ann[a], bjSeq, seq, flat );
					}
				}
			}
			bjSeqs.add( bjSeq );
		}
		return new SimpleSequenceIterator( bjSeqs.toArray( new org.biojava.bio.seq.Sequence[0] ) );
	}

	/**
	 * Returns a {@link RichSequence} for the {@link Sequence} <code>seq</code>
	 * trying to reconstruct accession, name, etc. from the
	 * {@link SequenceAnnotation} of <code>seq</code>.
	 * 
	 * @param bjAlph
	 *            the alphabet of the new {@link RichSequence}
	 * @param cont
	 *            the alphabet of <code>seq</code>
	 * @param seq
	 *            the {@link Sequence} to be converted
	 * @param defaultName
	 *            the default name, only used if no name could be found in the
	 *            annotation of <code>seq</code>
	 * 
	 * @return a new BioJava {@link RichSequence}
	 * 
	 * @throws BioException
	 *             if forwarded from BioJava
	 */
	private static RichSequence getSequence( FiniteAlphabet bjAlph, AlphabetContainer cont, Sequence seq, String defaultName ) throws BioException {
		SymbolList ssl = new SimpleSymbolList( bjAlph.getTokenization( "default" ), seq.toString() );
		Symbol[] symbols = (Symbol[])ssl.toList().toArray( new Symbol[0] );
		RichSequenceBuilder builder = (RichSequenceBuilder)RichSequenceBuilderFactory.FACTORY.makeSequenceBuilder();
		builder.startSequence();
		builder.addSymbols( bjAlph, symbols, 0, symbols.length );

		boolean setName = false;
		boolean setAccession = false;
		SequenceAnnotation[] ann = seq.getAnnotation();
		if( ann != null ) {
			int idx = -1;
			for( int i = 0; i < ann.length; i++ ) {
				if( ann[i].getIdentifier().equals( ANNOTATION_ID ) ) {
					idx = i;
					break;
				}
			}
			if( idx > -1 ) {

				Result[] res = ann[idx].getResults();

				for( int i = 0; i < res.length; i++ ) {
					if( res[i].getName().equals( "Description" ) ) {
						builder.setDescription( res[i].getResult().toString() );
					} else if( res[i].getName().equals( "Accession" ) ) {
						builder.setAccession( res[i].getResult().toString() );
						setAccession = true;
					} else if( res[i].getName().equals( "Name" ) ) {
						builder.setName( res[i].getResult().toString() );
						setName = true;
					} else if( res[i].getName().equals( "ID" ) ) {
						builder.setIdentifier( res[i].getResult().toString() );
					} else if( res[i].equals( "Version" ) ) {
						builder.setVersion( Integer.parseInt( res[i].getResult().toString() ) );
					}
				}
			}
		}
		if( !setName ) {
			builder.setName( defaultName );
		}
		if( !setAccession ) {
			builder.setAccession( defaultName );
		}
		builder.setNamespace( new SimpleNamespace( "default" ) );
		builder.endSequence();
		return builder.makeRichSequence();
	}

	/**
	 * Converts a {@link SequenceAnnotation} <code>annotation</code> and all of
	 * its sub-annotations to corresponding {@link Feature}s and adds the
	 * {@link Feature}s to <code>holder</code>.
	 * 
	 * @param annotation
	 *            the annotation of the sequence
	 * @param holder
	 *            the {@link FeatureHolder}, e.g. a
	 *            {@link org.biojava.bio.seq.Sequence} or a {@link Feature} *
	 * @param seq
	 *            the {@link Sequence}
	 * @param flat
	 *            indicates if the features should be flattened, this may be
	 *            necessary to preserve all features, because some data formats
	 *            do not support hierarchical features
	 * 
	 * @throws BioException
	 *             if forwarded from BioJava
	 * 
	 * @see SequenceAnnotation
	 */
	private static void annotationToFeatures( SequenceAnnotation annotation, FeatureHolder holder, Sequence seq, boolean flat ) throws BioException {
		Feature.Template templ = null;
		if( annotation instanceof StrandedLocatedSequenceAnnotationWithLength ) {
			templ = new StrandedFeature.Template();
			( (StrandedFeature.Template)templ ).strand = ( (StrandedLocatedSequenceAnnotationWithLength)annotation ).getStrandedness()
					.equals( StrandedLocatedSequenceAnnotationWithLength.Strand.FORWARD )	? StrandedFeature.POSITIVE
																							: StrandedFeature.NEGATIVE;
		}
		if( annotation instanceof LocatedSequenceAnnotationWithLength ) {
			templ = new Feature.Template();
			templ.location = new RangeLocation( ( (LocatedSequenceAnnotation)annotation ).getPosition() + 1,
					( (LocatedSequenceAnnotationWithLength)annotation ).getEnd() );
		} else if( annotation instanceof LocatedSequenceAnnotation ) {
			templ = new Feature.Template();
			templ.location = new PointLocation( ( (LocatedSequenceAnnotation)annotation ).getPosition() + 1 );
		} else {
			templ = new Feature.Template();
			templ.location = new RangeLocation( 1, seq.getLength() );
		}
		String t = annotation.getType();
		t = t.replace( ' ', '_' );
		if( t.length() > 15 ) {
			t = t.substring( 0, 15 );
		}
		templ.type = t;
		templ.source = "BioJavaAdapter: " + annotation.getType() + ", " + annotation.getIdentifier();
		Result[] res = annotation.getAnnotations();
		Map<String, String> map = new HashMap<String, String>();
		map.put( "identifier", annotation.getIdentifier() );
		if( res != null ) {
			for( int i = 0; i < res.length; i++ ) {
				map.put( res[i].getName(), res[i].getResult().toString() );
			}
		}
		templ.annotation = new SimpleAnnotation( map );
		Feature f = holder.createFeature( templ );
		SequenceAnnotation[] annotations = annotation.getSubAnnotations();
		if( annotations != null ) {
			for( int i = 0; i < annotations.length; i++ ) {
				if( flat ) {
					annotationToFeatures( annotations[i], holder, seq, flat );
				} else {
					annotationToFeatures( annotations[i], f, seq, flat );
				}
			}
		}
	}

	private static class WorkingWordTokenization extends NameTokenization {

		private String delim;

		private WorkingWordTokenization( FiniteAlphabet fab, boolean caseSensitive, String delim ) {
			super( fab, caseSensitive );
			this.delim = delim;
		}

		private WorkingWordTokenization( FiniteAlphabet fab, String delim ) {
			super( fab );
			this.delim = delim;
		}

		@Override
		public String tokenizeSymbolList( SymbolList sl ) throws IllegalSymbolException, IllegalAlphabetException {
			if( sl.getAlphabet() != getAlphabet() ) {
				throw new IllegalAlphabetException( "Alphabet " + sl.getAlphabet().getName() + " does not match " + getAlphabet().getName() );
			}
			StringBuffer sb = new StringBuffer();
			Iterator i = sl.iterator();
			while( i.hasNext() ) {
				Symbol sym = (Symbol)i.next();
				sb.append( tokenizeSymbol( sym ) );
				if( i.hasNext() ) {
					sb.append( delim );
				}
			}
			return sb.substring( 0 );
		}

		@Override
		protected List splitString( String seq ) throws IllegalSymbolException {
			String[] parts = seq.split( delim );
			ArrayList<String> list = new ArrayList<String>();
			for( int i = 0; i < parts.length; i++ ) {
				if( parts[i].length() > 0 ) {
					list.add( parts[i] );
				}
			}
			return list;
		}

	}
}

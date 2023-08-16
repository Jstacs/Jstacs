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

package de.jstacs.data.sequences.annotation;

import java.util.Collection;
import java.util.LinkedList;

import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.results.ListResult;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;

/**
 * Class for a general annotation of a {@link de.jstacs.data.sequences.Sequence}.
 * Annotations may be e.g. exons, introns, coding sequences or splice sites.
 * 
 * @author Jan Grau
 */
public class SequenceAnnotation extends ResultSet {

	private String type;

	private String identifier;

	/**
	 * Creates a new {@link SequenceAnnotation} of type <code>type</code> with
	 * identifier <code>identifier</code> and additional annotation (that does
	 * not fit the {@link SequenceAnnotation} definitions) given as a
	 * {@link Result} <code>result</code>.
	 * 
	 * @param type
	 *            the type of the annotation
	 * @param identifier
	 *            the identifier of the annotation
	 * @param result
	 *            the additional annotation
	 * 
	 * @see ResultSet#ResultSet(Result)
	 */
	public SequenceAnnotation( String type, String identifier, Result result ) {
		super( result );
		this.type = type;
		this.identifier = identifier;
	}

	/**
	 * Creates a new {@link SequenceAnnotation} of type <code>type</code> with
	 * identifier <code>identifier</code> and additional annotation (that does
	 * not fit the {@link SequenceAnnotation} definitions) given as an array of
	 * {@link Result}s <code>results</code>.
	 * 
	 * @param type
	 *            the type of the annotation
	 * @param identifier
	 *            the identifier of the annotation
	 * @param results
	 *            the additional annotation
	 * 
	 * @see ResultSet#ResultSet(Result[][])
	 */
	public SequenceAnnotation( String type, String identifier, Result[]... results ) {
		super( results );
		this.type = type;
		this.identifier = identifier;
	}

	/**
	 * Creates a new {@link SequenceAnnotation} of type <code>type</code> with
	 * identifier <code>identifier</code> and additional annotation (that does
	 * not fit the {@link SequenceAnnotation} definitions) given as an array of
	 * {@link Result}s <code>additionalAnnotation</code>. This
	 * {@link SequenceAnnotation} may contain sub-annotations
	 * <code>subAnnotations</code>, this may be e.g. the donor and acceptor site
	 * for splice sites or the exons for a gene.
	 * 
	 * @param type
	 *            the type of the annotation
	 * @param identifier
	 *            the identifier of the annotation
	 * @param subAnnotations
	 *            the sub-annotation
	 * @param additionalAnnotation
	 *            the additional annotation
	 * 
	 * @see ResultSet#ResultSet(Result[][])
	 */
	public SequenceAnnotation( String type, String identifier, SequenceAnnotation[] subAnnotations, Result... additionalAnnotation ) {
		super( fuse( additionalAnnotation, subAnnotations ) );
		this.type = type;
		this.identifier = identifier;
	}

	/**
	 * Creates a new {@link SequenceAnnotation} of type <code>type</code> with
	 * identifier <code>identifier</code> and additional annotation (that does
	 * not fit the {@link SequenceAnnotation} definitions) given as a
	 * {@link Collection} of {@link Result}s <code>results</code>.
	 * 
	 * @param type
	 *            the type of the annotation
	 * @param identifier
	 *            the identifier of the annotation
	 * @param results
	 *            the additional annotation
	 * 
	 * @see ResultSet#ResultSet(Collection)
	 */
	public SequenceAnnotation( String type, String identifier, Collection<? extends Result> results ) {
		super( results );
		this.type = type;
		this.identifier = identifier;
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link SequenceAnnotation} out of its XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link SequenceAnnotation} could not be reconstructed
	 *             out of the XML representation (the {@link StringBuffer}
	 *             <code>representation</code> could not be parsed)
	 * 
	 * @see ResultSet#ResultSet(StringBuffer)
	 * @see de.jstacs.Storable
	 */
	public SequenceAnnotation( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}

	/**
	 * Returns the type of this {@link SequenceAnnotation} as given in the
	 * constructor.
	 * 
	 * @return the type of this {@link SequenceAnnotation}
	 */
	public String getType() {
		return type;
	}

	/**
	 * Returns the additional annotations of this {@link SequenceAnnotation} as
	 * given in the constructor.
	 * 
	 * @return the additional annotations of this {@link SequenceAnnotation}
	 */
	public Result[] getAnnotations() {
		if( results.size() > 0 ) {
			LinkedList<Result> list = new LinkedList<Result>();
			for( int i = 0; i < results.size(); i++ ) {
				boolean add = true;
				if( results.get( i ) instanceof ListResult ) {
					ListResult lr = (ListResult)results.get( i );
					ResultSet[] rs = lr.getRawResult();
					for( int j = 0; j < rs.length; j++ ) {
						if( rs[j] instanceof SequenceAnnotation ) {
							add = false;
						}
					}
				}
				if( add ) {
					list.add( results.get( i ) );
				}
			}
			if( list.size() > 0 ) {
				return list.toArray( new Result[0] );
			} else {
				return null;
			}
		} else {
			return null;
		}
	}

	/**
	 * Returns the sub-annotations of this {@link SequenceAnnotation} as given
	 * in the constructor.
	 * 
	 * @return the sub-annotations of this {@link SequenceAnnotation}
	 */
	public SequenceAnnotation[] getSubAnnotations() {
		if( results.size() > 0 && results.get( results.size() - 1) instanceof ListResult ) {
			ListResult lr = (ListResult)results.get( results.size() - 1);
			ResultSet[] rs = lr.getRawResult();
			int num = 0;
			for( int i = 0; i < rs.length; i++ ) {
				if( rs[i] instanceof SequenceAnnotation ) {
					num++;
				}
			}
			if( num == 0 ) {
				return null;
			}
			SequenceAnnotation[] annot = new SequenceAnnotation[num];
			num = 0;
			for( int i = 0; i < rs.length; i++ ) {
				if( rs[i] instanceof SequenceAnnotation ) {
					annot[num] = (SequenceAnnotation)rs[i];
					num++;
				}
			}
			return annot;
		} else {
			return null;
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.results.ResultSet#toString()
	 */
	@Override
	public String toString() {
		StringBuffer buf = new StringBuffer();
		buf.append( type );
		buf.append( ", " );
		buf.append( identifier );
		buf.append( ":" );
		for( int i = 0; i < results.size(); i++ ) {
			if( results.get( i ) instanceof ListResult ) {
				ResultSet[] ress = (ResultSet[])results.get( i ).getValue();
				for( int j = 0; j < ress.length; j++ ) {
					buf.append( "\n" + ress[j].toString() );
				}
			} else {
				buf.append( "\n" + results.get( i ).toString() );
			}
		}
		return buf.toString();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.results.ResultSet#fromXML(java.lang.StringBuffer)
	 */
	@Override
	protected void fromXML( StringBuffer source ) throws NonParsableException {
		source = XMLParser.extractForTag( source, "sequenceAnnotation" );
		super.fromXML( XMLParser.extractForTag( source, "results" ) );
		type = XMLParser.extractObjectForTags( source, "type", String.class );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.results.ResultSet#toXML()
	 */
	@Override
	public StringBuffer toXML() {
		StringBuffer buf = super.toXML();
		XMLParser.addTags( buf, "results" );
		XMLParser.appendObjectWithTags( buf, type, "type" );
		XMLParser.addTags( buf, "sequenceAnnotation" );
		return buf;
	}

	/**
	 * Fuses the two arrays <code>res</code> and <code>res2</code> by writing
	 * the arrays one after the other into a new {@link Result} array, which is
	 * the super-type of {@link Result} as well as {@link SequenceAnnotation}.
	 * 
	 * @param res
	 *            the first array
	 * @param res2
	 *            the second array
	 * 
	 * @return the fused array
	 */
	private static Result[] fuse( Result[] res, SequenceAnnotation[] res2 ) {
		if( res2 == null || res2.length == 0 ) {
			return res;
		}
		Result[] ress = new Result[res.length + 1];
		for( int i = 0; i < res.length; i++ ) {
			ress[i] = res[i];
		}
		ress[ress.length - 1] = new ListResult( "Sub-annotation", "The list of sub-annotations", null, res2 );
		return ress;
	}

	/**
	 * Returns the identifier of this {@link SequenceAnnotation} as given in the
	 * constructor.
	 * 
	 * @return the identifier of this {@link SequenceAnnotation}
	 */
	public String getIdentifier() {
		return identifier;
	}

	@Override
	public int hashCode() {
		int hash = identifier.hashCode() + 31*type.hashCode();
		for(int i=0;i<results.size();i++){
			hash = 31*hash + results.get( i ).hashCode();
		}
		return hash;
	}
	
	

}

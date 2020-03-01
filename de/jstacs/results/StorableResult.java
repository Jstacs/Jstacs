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

package de.jstacs.results;

import de.jstacs.AnnotatedEntity;
import de.jstacs.DataType;
import de.jstacs.Storable;
import de.jstacs.classifiers.AbstractClassifier;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.SequenceScore;

/**
 * Class for {@link Result}s that are {@link Storable}s. The method
 * {@link #toXML()} is used to save the {@link StorableResult} together with the
 * result to an XML representation, or in the method {@link #toString()}.
 * 
 * @see de.jstacs.Storable
 * 
 * @author Jan Grau
 */
public class StorableResult extends Result {

	/**
	 * The {@link Storable} cannot be trained anyway.
	 */
	public static final byte NA = -1;
	/**
	 * The model/classifier has not been trained.
	 */
	public static final byte FALSE = 0;
	/**
	 * The model/classifier has been trained.
	 */
	public static final byte TRUE = 1;

	/**
	 * The {@link Storable} that is the result.
	 */
	private Storable object;

	/**
	 * Creates a result for an XML representation of an object.
	 * 
	 * @param name
	 *            the name of the result
	 * @param comment
	 *            a comment on the result
	 * @param object
	 *            the {@link Storable} that is the result
	 */
	public StorableResult( String name, String comment, Storable object ) {
		super( name, comment, DataType.STORABLE );
		this.object = object;
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a new {@link StorableResult} from its XML representation as
	 * returned by {@link #toXML()}.
	 * 
	 * @param buf
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if <code>buf</code> could not be parsed
	 */
	public StorableResult(StringBuffer buf) throws NonParsableException {
		super(buf);
	}

	/**
	 * Returns {@link #TRUE} if the model or classifier was trained when
	 * obtaining its XML representation stored in this {@link StorableResult},
	 * {@link #FALSE} if it was not, and {@link #NA} if the object could not be
	 * trained anyway.
	 * 
	 * @return if the model or classifier was trained or not
	 */
	public byte isInitialized() {
		if (object instanceof SequenceScore) {
			return ((SequenceScore) object).isInitialized() ? TRUE : FALSE;
		} else if (object instanceof AbstractClassifier) {
			return ((AbstractClassifier) object).isInitialized() ? TRUE : FALSE;
		} else {
			return NA;
		}
	}

	/**
	 * Returns the name of the class of the {@link Storable} corresponding to
	 * the XML representation stored in this {@link StorableResult}.
	 * 
	 * @return the name of the class
	 */
	public String getClassName() {
		return object.getClass().getName();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.results.Result#getResult()
	 */
	@Override
	public String getValue() {
		return object.toXML().toString();
	}

	/**
	 * Returns the instance of the {@link Storable} that is the result of this
	 * {@link StorableResult}.
	 * 
	 * @return the instance
	 */
	public Storable getResultInstance() {
		return object;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.AnnotatedEntity#getXMLTag()
	 */
	@Override
	public String getXMLTag() {
		return "objectResult";
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.AnnotatedEntity#appendFurtherInfos(java.lang.StringBuffer)
	 */
	@Override
	protected void appendFurtherInfos( StringBuffer buf ) {
		XMLParser.appendObjectWithTags( buf, object, "object" );
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.AnnotatedEntity#extractFurtherInfos(java.lang.StringBuffer)
	 */
	@Override
	protected void extractFurtherInfos( StringBuffer representation ) throws NonParsableException {
		this.object = XMLParser.extractObjectForTags( representation, "object", Storable.class );
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return name + ":\n" + object.toXML() + "\n";
	}

	@Override
	public boolean equalValues(AnnotatedEntity a) {
		if(a instanceof StorableResult) {
			StorableResult s = (StorableResult)a;
			if(s == this) {
				return true;
			}else if(s.object == this.object) {
				return true;
			}else if(s.object != null && this.object != null && s.object.equals(this.object)) {
				return true;
			}else {
				return false;
			}
		}else {
			return false;
		}
	}
	
	
	
	
}

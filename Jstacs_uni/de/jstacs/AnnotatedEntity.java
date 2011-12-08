package de.jstacs;

import de.jstacs.io.XMLParser;


/**
 * Superclass for all Jstacs entities that have a name, a comment, and a data type as annotations.
 * Specifically, such entities are {@link de.jstacs.results.Result}s for storing results of some computation together with an annotation on their
 * meaning, and {@link de.jstacs.parameters.Parameter}s, which allow for annotating external parameters, especially parameters of constructors.
 * 
 * @author Jan Grau, Jens Keilwagen
 *
 */
public abstract class AnnotatedEntity implements Storable {

	/**
	 * The name of the entity.
	 */
	protected String name;

	/**
	 * The comment for the entity.
	 */
	protected String comment;

	/**
	 * The data type of the entity.
	 */
	protected DataType datatype;

	/**
	 * The main constructor which takes the main information of a {@link AnnotatedEntity}.
	 * 
	 * @param name
	 *            the name of the result
	 * @param comment
	 *            the comment for the result
	 * @param datatype
	 *            the data type of the result
	 */
	protected AnnotatedEntity(String name, String comment, DataType datatype) {
		this.name = name;
		this.comment = comment;
		this.datatype = datatype;
	}

	/**
	 * The standard constructor for the interface {@link Storable}. Creates a
	 * new {@link AnnotatedEntity} out of its XML representation.
	 * 
	 * @param rep
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML representation is not parsable
	 * 
	 * @see Storable
	 * @see #extractFurtherInfos(StringBuffer)
	 */
	protected AnnotatedEntity(StringBuffer rep) throws NonParsableException {
		rep = XMLParser.extractForTag( rep, getXMLTag() );
		name = XMLParser.extractObjectForTags( rep, "name", String.class );
		comment = XMLParser.extractObjectForTags( rep, "comment", String.class );
		datatype = XMLParser.extractObjectForTags( rep, "datatype", DataType.class );
		extractFurtherInfos( rep );
	}
	
	/**
	 * This method returns a tag used as outer tag of the XML description.
	 * @return a tag used as outer tag of the XML description
	 */
	public abstract String getXMLTag();
		
	public final StringBuffer toXML() {
		StringBuffer buf = new StringBuffer();
		XMLParser.appendObjectWithTags( buf, name, "name" );
		XMLParser.appendObjectWithTags( buf, comment, "comment" );
		XMLParser.appendObjectWithTags( buf, datatype, "datatype" );
		appendFurtherInfos( buf );
		return buf;
	}
	/**
	 * This method can be used in the method {@link Storable#toXML()} to extract
	 * further information (name, comment, datatype).
	 * 
	 * @param buf
	 *            a XML representation of the main information as
	 *            {@link StringBuffer}
	 * 
	 * @see Storable#toXML()
	 */
	protected abstract void appendFurtherInfos( StringBuffer buf );

	/**
	 * This method can be used in the constructor with parameter {@link StringBuffer} to
	 * extract the further information.
	 * 
	 * @param buf
	 *            a XML represenation of the main information as
	 *            {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML representation is not parsable
	 * 
	 * @see #AnnotatedEntity(StringBuffer)
	 */
	protected abstract void extractFurtherInfos(StringBuffer buf) throws NonParsableException;
	
	
	/**
	 * Returns the data type of the {@link AnnotatedEntity}.
	 * 
	 * @return the data type of the {@link AnnotatedEntity}
	 */
	public final DataType getDatatype() {
		return datatype;
	}

	/**
	 * Returns the value of the {@link AnnotatedEntity}.
	 * 
	 * @return the value of the {@link AnnotatedEntity}
	 */
	public abstract Object getValue();

	/**
	 * Returns the name of the {@link AnnotatedEntity}.
	 * 
	 * @return the name of the {@link AnnotatedEntity}
	 */
	public final String getName() {
		return name;
	}

	/**
	 * Returns the comment on the {@link AnnotatedEntity}.
	 * 
	 * @return the comment on the {@link AnnotatedEntity}
	 */
	public final String getComment() {
		return comment;
	}
}
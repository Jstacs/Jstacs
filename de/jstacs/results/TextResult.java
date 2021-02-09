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

package de.jstacs.results;

import de.jstacs.DataType;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.FileParameter.FileRepresentation;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.results.savers.TextResultSaver;
import de.jstacs.tools.ui.galaxy.Galaxy;


/**
 * Class for a result that is basically a text file (or its contents).
 * This is the {@link Result} analogon of a {@link FileParameter}, but conceptually (though not formally)
 * restricted to text contents.
 * 
 * @author Jan Grau
 *
 */
public class TextResult extends Result {

	private FileRepresentation value;
	private String mime;
	private String extendedType;
	private String producer;
	private boolean export;
	private boolean isTempFile;
	
	/**
	 * Creates a new {@link TextResult} with given name, comment, content, mime type, and additional info, assuming
	 * that <code>file</code> is a temporary file (if it has a filename associated and this file exists).
	 * @param name the name of the result
	 * @param comment a comment on the result
	 * @param file the contents of the result, encapsulated in a {@link FileRepresentation}.
	 * @param mime the mime type of the content
	 * @param producer the producer (may be <code>null</code>)
	 * @param extendedType the extended type (may be <code>null</code>)
	 * @param export if <code>true</code>, the contents are exported as a separate file when used in {@link Galaxy}
	 */
	public TextResult(String name, String comment, FileRepresentation file, String mime, String producer, String extendedType, boolean export){
		super(name, comment, DataType.FILE);
		this.value = file;
		this.mime = mime;
		this.producer = producer;
		this.extendedType = extendedType;
		this.export = export;
		this.isTempFile = true;
	}
	
	/**
	 * Creates a new {@link TextResult} with given name, comment, content, mime type, and additional info.
	 * @param name the name of the result
	 * @param comment a comment on the result
	 * @param file the contents of the result, encapsulated in a {@link FileRepresentation}.
	 * @param isTempFile if <code>file</code> refers to a temporary file. If <code>false</code>, the corresponding original file will not be removed by {@link TextResultSaver#writeOutput(TextResult, java.io.File)} 
	 * @param mime the mime type of the content
	 * @param producer the producer (may be <code>null</code>)
	 * @param extendedType the extended type (may be <code>null</code>)
	 * @param export if <code>true</code>, the contents are exported as a separate file when used in {@link Galaxy}
	 */
	public TextResult(String name, String comment, FileRepresentation file, boolean isTempFile, String mime, String producer, String extendedType, boolean export){
		super(name, comment, DataType.FILE);
		this.value = file;
		this.mime = mime;
		this.producer = producer;
		this.extendedType = extendedType;
		this.export = export;
		this.isTempFile = isTempFile;
	}
	
	/**
	 * Creates a {@link TextResult} from its XML representation.
	 * @param xml the XML representation
	 * @throws NonParsableException if XML could not be parsed
	 */
	public TextResult(StringBuffer xml) throws NonParsableException{
		super(xml);
	}
	
	/**
	 * Returns <code>true</code> if the contents are saved to a separate file in {@link Galaxy}.
	 * @return if the contents are saved to a separate file
	 */
	public boolean getExport(){
		return export;
	}
	
	/**
	 * Checks if the list of mime types given in <code>p1</code>
	 * contains an element that is equal to one of the mime types given in <code>mime2</code> (may be multiple, separated by commas).
	 * @param p1 the reference mime types
	 * @param mime2 the mime types to be checked
	 * @return if <code>p1</code> and <code>mime2</code> share a mime type
	 */
	public static boolean equals(String[] p1, String mime2){
		if(p1 == null || mime2 == null){
			return false;
		}
		String[] p2 = mime2.split( "\\," );
		for(int i=0;i<p1.length;i++){
			for(int j=0;j<p2.length;j++){
				if(p1[i].equals( p2[j] )){
					return true;
				}
			}
		}
		return false;
	}
	
	/**
	 * Checks if the list of mime types given in <code>mime1</code> (may be multiple, separated by commas)
	 * contains an element that is equal to one of the mime types given in <code>mime2</code>.
	 * @param mime1 the reference mime types
	 * @param mime2 the mime types to be checked
	 * @return if <code>mime1</code> and <code>mime2</code> share a mime type
	 */
	public static boolean equals(String mime1, String mime2){
		if(mime1 == null || mime2 == null){
			return false;
		}
		String[] p1 = mime1.split( "\\," );
		return equals( p1, mime2 );
	}
	
	/**
	 * Fills the supplied {@link FileParameter} with a clone of the contents of this {@link TextResult}.
	 * @param par the {@link FileParameter} to be filled
	 * @throws IllegalValueException if the value of this {@link FileParameter} is not accepted by the supplied parameter
	 * @throws CloneNotSupportedException if the contents of this {@link TextResult} could not be cloned.
	 */
	public void fill(FileParameter par) throws IllegalValueException, CloneNotSupportedException {
		if(equals( par.getAcceptedMimeType(), mime )){
			FileRepresentation temp = value.clone();
			if(temp.getFilename() == null || temp.getFilename().length() == 0){
				temp.setFilename(this.getName());
			}
			if(temp.getExtension() == null){
				temp.setExtension( this.mime );
			}
			par.setValue( temp );
		}else{
			throw new IllegalValueException(name,"Mime "+mime+" does not match required mime "+par.getAcceptedMimeType()+".");
		}
	}
	
	/**
	 * Returns the mime type of this {@link TextResult}.
	 * @return the mime type
	 */
	public String getMime() {
		return mime;
	}
	
	/**
	 * Returns the extended type of this {@link TextResult}. May be used for checking compatibility for
	 * {@link TextResult#fill(FileParameter)} using {@link FileParameter#getExtendedType()}.
	 * @return the extended type
	 */
	public String getExtendedType() {
		return extendedType;
	}

	
	/**
	 * Sets the extended type of this {@link TextResult}.
	 * @param extendedType the extended type
	 * @see FileParameter#setExtendedType(String)
	 */
	public void setExtendedType( String extendedType ) {
		this.extendedType = extendedType;
	}

	public boolean isTempFile() {
		return isTempFile;
	}
	
	public void setTempFile(boolean set) {
		this.isTempFile = set;
	}
	
	@Override
	public String getXMLTag() {
		return "TextResult";
	}

	@Override
	protected void appendFurtherInfos( StringBuffer buf ) {
		XMLParser.appendObjectWithTags(buf, mime, "mime");
		if(isTempFile()) {
			value.getContent();
			value.setFilename("");
		}
		XMLParser.appendObjectWithTags(buf, value,"value");
		XMLParser.appendObjectWithTags(buf, producer,"producer");
		XMLParser.appendObjectWithTags(buf, export, "export");
		XMLParser.appendObjectWithTags( buf, extendedType, "extype" );
		XMLParser.appendObjectWithTags(buf, isTempFile, "isTempFile");
	}

	@Override
	protected void extractFurtherInfos( StringBuffer buf ) throws NonParsableException {
		mime = (String)XMLParser.extractObjectForTags( buf, "mime" );
		value = (FileRepresentation)XMLParser.extractObjectForTags( buf, "value" );
		producer = (String)XMLParser.extractObjectForTags( buf, "producer" );
		try{
			export = (Boolean) XMLParser.extractObjectForTags(buf, "export");
		}catch(NonParsableException e){
			export = false;
		}
		extendedType = (String)XMLParser.extractObjectForTags( buf, "extype" );
		try {
			isTempFile = (Boolean) XMLParser.extractObjectForTags(buf, "isTempFile");
		}catch(NonParsableException e) {
			isTempFile = true;
		}
	}

	@Override
	public FileRepresentation getValue() {
		return value;
	}


	/**
	 * Returns the producer (i.e., the tool/application/program) that created this {@link TextResult}.
	 * @return the producer
	 */
	public Object getProducer() {
		return producer;
	}

}

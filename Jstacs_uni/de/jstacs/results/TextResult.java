package de.jstacs.results;

import de.jstacs.DataType;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.FileParameter.FileRepresentation;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;


public class TextResult extends Result {

	private FileRepresentation value;
	private String mime;
	private String extendedType;
	private String producer;
	private boolean export;
	
	public TextResult(String name, String comment, FileRepresentation file, String mime, String producer, String extendedType, boolean export){
		super(name, comment, DataType.FILE);
		this.value = file;
		this.mime = mime;
		this.producer = producer;
		this.extendedType = extendedType;
		this.export = export;
	}
	
	public TextResult(StringBuffer xml) throws NonParsableException{
		super(xml);
	}
	
	public boolean getExport(){
		return export;
	}
	
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
	
	public static boolean equals(String mime1, String mime2){
		if(mime1 == null || mime2 == null){
			return false;
		}
		String[] p1 = mime1.split( "\\," );
		return equals( p1, mime2 );
	}
	
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
			throw new IllegalValueException("Mime "+mime+" does not match required mime "+par.getAcceptedMimeType()+".");
		}
	}
	
	public String getMime() {
		return mime;
	}
	
	
	public String getExtendedType() {
		return extendedType;
	}

	
	public void setExtendedType( String extendedType ) {
		this.extendedType = extendedType;
	}

	@Override
	public String getXMLTag() {
		return "TextResult";
	}

	@Override
	protected void appendFurtherInfos( StringBuffer buf ) {
		XMLParser.appendObjectWithTags(buf, mime, "mime");
		XMLParser.appendObjectWithTags(buf, value,"value");
		XMLParser.appendObjectWithTags(buf, producer,"producer");
		XMLParser.appendObjectWithTags( buf, extendedType, "extype" );
	}

	@Override
	protected void extractFurtherInfos( StringBuffer buf ) throws NonParsableException {
		mime = (String)XMLParser.extractObjectForTags( buf, "mime" );
		value = (FileRepresentation)XMLParser.extractObjectForTags( buf, "value" );
		producer = (String)XMLParser.extractObjectForTags( buf, "producer" );
		extendedType = (String)XMLParser.extractObjectForTags( buf, "extype" );
	}

	@Override
	public FileRepresentation getValue() {
		return value;
	}


	public Object getProducer() {
		return producer;
	}

}

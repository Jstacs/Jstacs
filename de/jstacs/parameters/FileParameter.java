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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.parameters;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.StringReader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;

import de.jstacs.DataType;
import de.jstacs.Storable;
import de.jstacs.io.FileManager;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.parameters.validation.ParameterValidator;
import de.jstacs.tools.ui.cli.CLI;
import de.jstacs.tools.ui.galaxy.GalaxyAdaptor;
import de.jstacs.utils.Compression;

/**
 * Class for a {@link Parameter} that represents a local file.
 * 
 * @author Jan Grau
 * 
 */
public class FileParameter extends Parameter implements GalaxyConvertible {

	
	/**
	 * <code>true</code> if the parameter is required
	 */
	private boolean required;
	/**
	 * The MIME-type of accepted files
	 */
	private String mime;
	
	/**
	 * More specific file type of accepted files
	 * for programmatic use
	 */
	private String extendedType;
	/**
	 * The file
	 */
	private FileRepresentation value;
	/**
	 * The default value
	 */
	private FileRepresentation defaultValue;
	/**
	 * <code>true</code> if a file is set as value
	 */
	private boolean isSet;
	/**
	 * The error message, <code>null</code> if no error occurred
	 */
	private String errorMessage;

	/**
	 * The parameter validator
	 */
	private ParameterValidator valid;

	/**
	 * if true the MIME type is checked in {@link CLI}
	 */
	private boolean checkMimeType;
	
	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#clone()
	 */
	@Override
	public FileParameter clone() throws CloneNotSupportedException {
		FileParameter clone = (FileParameter) super.clone();
		clone.value = value == null ? null : value.clone();
		clone.defaultValue = defaultValue == null ? null : defaultValue.clone();
		return clone;
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Restores a {@link FileParameter} from an XML representation.
	 * 
	 * @param buf
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML code could not be parsed
	 */
	public FileParameter(StringBuffer buf) throws NonParsableException {
		super( buf );
	}

	/**
	 * Creates a {@link FileParameter}.
	 * 
	 * @param name
	 *            the name of the parameter
	 * @param comment
	 *            a comment on the parameter
	 * @param filetype
	 *            the type of allowed files, may be allowed file extensions, separated by commas
	 * @param required
	 *            <code>true</code> if this {@link FileParameter} is required to
	 *            continue, <code>false</code> otherwise
	 */
	public FileParameter(String name, String comment, String filetype, boolean required) {
		super( name, comment, DataType.FILE );
		this.mime = filetype;
		this.required = required;
	}

	/**
	 * Constructs a {@link FileParameter}.
	 * 
	 * @param name
	 *            the name of the parameter
	 * @param comment
	 *            a comment on the parameter
	* @param filetype
	 *            the type of allowed files, may be allowed file extensions, separated by commas
	 * @param required
	 *            <code>true</code> if this {@link FileParameter} is required
	 * @param validator
	 *            a validator that validates e.g. the contents of the file
	 */
	public FileParameter(String name, String comment, String filetype, boolean required, ParameterValidator validator) {
		this(name, comment, filetype, required, validator, false);
	}
	
	/**
	 * Constructs a {@link FileParameter}.
	 * 
	 * @param name
	 *            the name of the parameter
	 * @param comment
	 *            a comment on the parameter
	 * @param filetype
	 *            the type of allowed files, may be allowed file extensions, separated by commas
	 * @param required
	 *            <code>true</code> if this {@link FileParameter} is required
	 * @param validator
	 *            a validator that validates e.g. the contents of the file
	 * @param checkMimeType
	 *            allowing to switch between checking mime type or not if values are set using {@link CLI}
	 */
	public FileParameter(String name, String comment, String filetype, boolean required, ParameterValidator validator, boolean checkMimeType) {
		this(name, comment, filetype, required);
		this.valid = validator;
		this.checkMimeType=checkMimeType;
	}

	/**
	 * Returns <code>true</code> if the MIME type should be checked
	 * 
	 * @return <code>true</code> if the MIME type should be checked
	 * 
	 * @see CLI
	 */
	public boolean checkMimeType() {
		return checkMimeType;
	}
	
	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#isAtomic()
	 */
	@Override
	public boolean isAtomic() {
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#isRequired()
	 */
	@Override
	public boolean isRequired() {
		return required;
	}

	/**
	 * Resets the {@link FileParameter} to its original state.
	 */
	@Override
	public void reset() {
		this.value = null;
		this.isSet = false;
		this.errorMessage = null;
	}

	/**
	 * Returns the content of the file.
	 * 
	 * @return the content of the file
	 */
	public FileRepresentation getFileContents() {
		return value;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#getErrorMessage()
	 */
	@Override
	public String getErrorMessage() {
		return errorMessage;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#checkValue(java.lang.Object)
	 */
	@Override
	public boolean checkValue(Object value) {
		if(value instanceof String){
			value = new FileRepresentation( (String)value );
		}
		if (valid != null) {
			
			if (valid.checkValue(value)) {
				errorMessage = null;
				return true;
			} else {
				errorMessage = valid.getErrorMessage();
				return false;
			}
		} else if (value != null && value instanceof FileRepresentation) {

			FileRepresentation f = (FileRepresentation) value;
			if (f.getFilename() != null && f.getFilename().length() != 0
					|| f.getContent() != null && f.getContent().length() != 0) {
				errorMessage = null;
				return true;
			} else {
				errorMessage = "No file specified or file is empty.";
				return false;
			}
		} else {
			errorMessage = "Value is no file or null.";
			return false;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#setDefault(java.lang.Object)
	 */
	@Override
	public void setDefault(Object defaultValue) throws IllegalValueException {
		if(defaultValue instanceof String){
			defaultValue = new FileRepresentation( (String)defaultValue );
		}
		if (checkValue(defaultValue)) {
			this.defaultValue = (FileRepresentation) defaultValue;
			setValue(defaultValue);
			isSet = false;
		} else {
			throw new IllegalValueException(name,errorMessage);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#setValue(java.lang.Object)
	 */
	@Override
	public void setValue(Object value) throws IllegalValueException {
		if(value instanceof String){
			value = new FileRepresentation( (String)value );
		}
		if (!checkValue(value)) {
			value = null;
			isSet = false;
			throw new IllegalValueException(name,errorMessage);
		}
		this.value = (FileRepresentation) value;
		this.isSet = true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#getValue()
	 */
	@Override
	public String getValue() {
		if (value == null) {
			return null;
		} else {
			return value.getFilename();
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#hasDefaultOrIsSet()
	 */
	@Override
	public boolean hasDefaultOrIsSet() {
		return isSet();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#isSet()
	 */
	@Override
	public boolean isSet() {
		return isSet;
	}

	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.AnnotatedEntity#getXMLTag()
	 */
	@Override
	public String getXMLTag() {
		return "fileParameter";
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.parameters.Parameter#appendFurtherInfos(java.lang.StringBuffer)
	 */
	@Override
	protected void appendFurtherInfos( StringBuffer buf ) {
		super.appendFurtherInfos( buf );
		
		XMLParser.appendObjectWithTags(buf, mime, "mime");
		XMLParser.appendObjectWithTags(buf, checkMimeType, "checkMimeType");
		XMLParser.appendObjectWithTags(buf, required, "required");
		XMLParser.appendObjectWithTags(buf, isSet, "isSet");
		XMLParser.appendObjectWithTags(buf, errorMessage, "errorMessage");
		
		XMLParser.appendObjectWithTags(buf, value,"value");
		XMLParser.appendObjectWithTags(buf, valid, "validator");
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.parameters.Parameter#extractFurtherInfos(java.lang.StringBuffer)
	 */
	@Override
	protected void extractFurtherInfos( StringBuffer buf ) throws NonParsableException {
		super.extractFurtherInfos( buf );
		
		mime = XMLParser.extractObjectForTags(buf, "mime", String.class );
		checkMimeType = XMLParser.extractObjectForTags(buf, "checkMimeType", boolean.class );
		required = XMLParser.extractObjectForTags(buf, "required", boolean.class );
		isSet = XMLParser.extractObjectForTags(buf, "isSet", boolean.class );
		errorMessage = XMLParser.parseString( XMLParser.extractObjectForTags(buf, "errorMessage", String.class ) );
		
		value = XMLParser.extractObjectForTags(buf, "value", FileRepresentation.class );
		valid = XMLParser.extractObjectForTags(buf, "validator", ParameterValidator.class );
	}

	/**
	 * Returns the type(s) of the allowed files as list of file extensions, separated by commas
	 * 
	 * @return the type(s) of the allowed files
	 */
	public String getAcceptedMimeType() {
		return mime;
	}
	
	/**
	 * Returns the extended type (or <code>null</code> if not set) of this {@link FileParameter}.
	 * @return the extended type
	 */
	public String getExtendedType() {
		return extendedType;
	}

	/**
	 * Sets the extended type of this {@link FileParameter}.
	 * @param extendedType the extended type
	 */
	public void setExtendedType( String extendedType ) {
		this.extendedType = extendedType;
	}

	private String mimeToGalaxy(String mime){
		String[] supported = new String[]{"gb","genbank","ab1", "afg", "arff", "asn1", "asn1-binary", "axt", "fli", "bam", "bed", "bedgraph", "bedstrict", "bed6", "bed12", "len", "bigbed", "bigwig", "cxb", "chrint", "csv", "customtrack", "bowtie_color_index", "bowtie_base_index", "csfasta", "data", "data_manager_json", "fasta", "fasta.gz", "fastq", "fastq.gz", "fastqsanger", "fastqsolexa", "fastqcssanger", "fastqillumina", "fqtoc", "eland", "elandmulti", "genetrack", "gff", "gff3", "gif", "gmaj.zip", "gtf", "toolshed.gz", "h5", "html", "interval", "picard_interval_list", "gatk_interval", "gatk_report", "gatk_dbsnp", "gatk_tranche", "gatk_recal", "jpg", "tiff", "bmp", "im", "pcd", "pcx", "ppm", "psd", "xbm", "xpm", "rgb", "pbm", "pgm", "rna_eps", "searchgui_archive", "peptideshaker_archive", "eps", "rast", "laj", "lav", "maf", "mafcustomtrack", "encodepeak", "pdf", "pileup", "obo", "owl", "png", "qual", "qualsolexa", "qualillumina", "qualsolid", "qual454", "Roadmaps", "sam", "scf", "Sequences", "snpeffdb", "snpsiftdbnsfp", "dbnsfp.tabular", "sff", "svg", "taxonomy", "tabular", "twobit", "sqlite", "gemini.sqlite", "txt", "linecount", "memexml", "cisml", "xml", "vcf", "bcf", "velvet", "wig", "interval_index", "tabix", "bgzip", "vcf_bgzip", "bgzip", "phyloxml", "nhx", "nex", "affybatch", "eigenstratgeno", "eigenstratpca", "eset", "fped", "fphe", "gg", "ldindep", "malist", "lped", "pbed", "pheno", "pphe", "rexpbase", "rgenetics", "snptest", "snpmatrix", "xls", "ipynb", "json", "xgmml", "sif", "rdf", "xlsx"};
		HashSet<String> set = new HashSet<String>();
		Collections.addAll( set, supported );
		String[] mimes = mime.split( "," );
		LinkedList<String> res = new LinkedList<String>();
		for(int i=0;i<mimes.length;i++){
			mimes[i] = mimes[i].trim();
			if(set.contains( mimes[i] )){
				res.add( mimes[i] );
			}
		}
		if(res.size() == 0){
			return mime;
		}else{
			StringBuffer sb = new StringBuffer();
			for(int i=0;i<res.size();i++){
				sb.append( res.get( i ) );
				if(i<res.size()-1){
					sb.append( "," );
				}
			}
			return sb.toString();
		}
	}
	
	@Override
	public void toGalaxy( String namePrefix, String configPrefix, int depth, StringBuffer descBuffer, StringBuffer configBuffer, boolean addLine, int indentation  ) {
		namePrefix = namePrefix+"_"+GalaxyAdaptor.getLegalName( getName() );
		StringBuffer buf = new StringBuffer();
		String line = "";
		if(addLine){
			//line = "&lt;hr style=&quot;height:2px;background-color:"+GalaxyAdaptor.getColor( depth )+";color:"+GalaxyAdaptor.getColor( depth )+";border:none&quot; /&gt;";
			line = "&lt;hr /&gt;";
		}
		XMLParser.addTagsAndAttributes( buf, "param", "type=\"data\" format=\""+mimeToGalaxy( mime )+"\" name=\""+namePrefix+"\" label=\""+line+getName()+"\" help=\""+getComment()+"\" value=\""+(defaultValue == null ? "" : defaultValue)+"\" optional=\""+(!isRequired())+"\"", indentation );
		descBuffer.append( buf );
		buf = new StringBuffer();
		buf.append( "${"+configPrefix+namePrefix+"}" );
		XMLParser.addTags( buf, "value" );
		StringBuffer buf2 = new StringBuffer();
		buf2.append( "${"+configPrefix+namePrefix+".ext}" );
		XMLParser.addTags( buf2, "extension" );
		buf.append( buf2 );
		XMLParser.addTags( buf, namePrefix );
		configBuffer.append( buf );
	}

	public String toString(){
		return name + " (" + comment
				+ (defaultValue!=null?", default = " + defaultValue.getFilename():"")
				+ (mime!=null?", type = " + mime:"")
				+ (required ? "" : ", OPTIONAL" )
				+ ")\t= " + (value != null ? value.getFilename() : "null");
	}
	
	@Override
	public void fromGalaxy( String namePrefix, StringBuffer command ) throws Exception {
		namePrefix = namePrefix+"_"+GalaxyAdaptor.getLegalName( getName() );
		StringBuffer cont = XMLParser.extractForTag( command, namePrefix );
		try{
			String val = XMLParser.extractForTag( cont, "value" ).toString();
			if( val.equalsIgnoreCase("none") ) {
				this.isSet=false;
				return;
			}
			String ext = XMLParser.extractForTag( cont, "extension" ).toString().trim();
			this.value = new FileRepresentation( val );
			this.value.setExtension( ext );
		}catch(NonParsableException e){
			this.value = new FileRepresentation( cont.toString() );
		}
		this.isSet = true;
	}
	
	@Override
	public void toGalaxyTest(String namePrefix, int depth, StringBuffer testBuffer, int indentation) throws Exception {
		if( isSet ) {
			XMLParser.insertIndentation(testBuffer, indentation, testBuffer.length());
			testBuffer.append("<param name=\"" + namePrefix+"_"+GalaxyAdaptor.getLegalName( getName() ) + "\" file=\"" + (String)getValue() + "\" />\n");
		}
	}
	
	/**
	 * Class that represents a file.
	 * 
	 * @author Jan Grau, Jens Keilwagen
	 */
	public static class FileRepresentation implements Storable, Cloneable {

		/**
		 * The name of the file
		 */
		private String filename;
		
		private String ext;
		
		/**
		 * The contents of the file
		 */
		private String content;

		private boolean compressed;

		/**
		 * Creates a {@link FileRepresentation} out of the filename and the
		 * file's contents.
		 * 
		 * @param filename
		 *            the name of the file
		 * @param content
		 *            the contents of the file
		 */
		public FileRepresentation(String filename, String content) {
			this.filename = filename;
			
			if(content.length() < 100000){
				this.content = content;
			}else{
				try{
					this.content = Compression.zip( content );
					this.compressed = true;
				}catch(IOException e){
					this.content = content;
					this.compressed = false;
				}
			}
			//this.content = content;
			int idx = filename.lastIndexOf( '.' );
			if(idx >= 0){
				this.ext = filename.substring( idx+1 );
			}
		}
		
		/**
		 * Creates a new {@link FileRepresentation} from a filename.
		 * @param filename the filename
		 */
		public FileRepresentation(String filename){
			this.filename = filename;
			int idx = filename.lastIndexOf( '.' );
			if(idx >= 0){
				this.ext = filename.substring( idx+1 );
			}
		}

		/**
		 * The standard constructor for the interface {@link de.jstacs.Storable}
		 * . Restores the {@link FileRepresentation} from an XML representation.
		 * 
		 * @param buf
		 *            the XML representation as {@link StringBuffer}
		 * 
		 * @throws NonParsableException
		 *             if the {@link StringBuffer} could not be parsed
		 */
		public FileRepresentation(StringBuffer buf) throws NonParsableException {
			fromXML(buf);
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Object#clone()
		 */
		@Override
		public FileRepresentation clone() throws CloneNotSupportedException {
			return (FileRepresentation) super.clone();
		}

		/**
		 * Returns the filename.
		 * 
		 * @return the name of the file
		 */
		public String getFilename() {
			return filename;
		}
		
		/**
		 * Sets the file name of this {@link FileParameter}
		 * @param filename the new file name
		 */
		public void setFilename(String filename) {
			this.filename = filename;
		}
		
		/**
		 * Returns the size of the file of the {@link FileRepresentation}
		 * as specified in the constructor. If  <code>filename</code> is <code>null</code> in
		 * the constructor, -1 is returned.
		 * @return the size of the file
		 * @see File#length()
		 */
		public long getFilesize(){
			if(filename == null){
				return -1;
			}else{
				return (new File(filename)).length();
			}
		}

		/**
		 * Returns the content of the file.
		 * 
		 * @return the content of the file
		 */
		public String getContent() {
			if(content == null){
				if(filename != null){
					try{
						String temp = FileManager.readFile(filename).toString();
						if(temp.length() < 100000){
							this.content = temp;
						}else{
							try{
								this.content = Compression.zip( temp );
								this.compressed = true;
							}catch(IOException e){
								this.content = temp;
								this.compressed = false;
							}
						}
						return temp;
					}catch(IOException e){
						return null;
					}
				}
			}else if(compressed){
				try {
					return Compression.unzip( content );
				} catch ( Exception e ) {
					e.printStackTrace();
					return null;
				}
			}
			return content;
		}
		
		/**
		 * Sets the extension of this {@link FileParameter}.
		 * @param ext the extension
		 */
		public void setExtension(String ext){
			this.ext = ext;
		}
		
		/**
		 * Returns the extension of this {@link FileParameter}.
		 * @return the extension
		 */
		public String getExtension(){
			return ext;
		}
		
		/**
		 * This method moves the underlying file to the specified location. If the underlying file does not exist, it return false.
		 * 
		 * @param newFileName the new file name
		 * 
		 * @return <code>true</code>
		 * 
		 * @throws IOException if something went wrong while moving the file
		 * 
		 * @see Files#move(java.nio.file.Path, java.nio.file.Path, java.nio.file.CopyOption...)
		 */
		public boolean moveFile( String newFileName ) throws IOException {
			boolean moved = false;
			if( filename != null ) {
				Files.move( Paths.get(filename), Paths.get(newFileName), StandardCopyOption.REPLACE_EXISTING );
				setFilename(newFileName);
				moved=true;
			}
			return moved;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see de.jstacs.Storable#toXML()
		 */
		public StringBuffer toXML() {
			StringBuffer buf = new StringBuffer();
			if( filename == null || !new File( filename ).exists() ) {
				if( content == null ) {
					System.err.println("WARNING: empty TextResult"); //TODO
				}
				filename=null;
				XMLParser.appendObjectWithTags(buf, content, "content");
				XMLParser.appendObjectWithTags(buf, ext, "ext");
			}
			XMLParser.appendObjectWithTags(buf, filename, "filename");
			XMLParser.appendObjectWithTags( buf, compressed, "compressed" );
			XMLParser.addTags(buf, "fileRepresentation");

			return buf;
		}

		private void fromXML(StringBuffer representation)
				throws NonParsableException {
			representation = XMLParser.extractForTag(representation,
					"fileRepresentation");
			filename = XMLParser.extractObjectForTags(representation, "filename", String.class );
			compressed = XMLParser.extractObjectForTags( representation, "compressed", Boolean.class );
			if(filename != null) {
				int idx = filename.lastIndexOf( '.' );
				if(idx >= 0){
					this.ext = filename.substring( idx+1 );
				}
			} else {
				try{ 
					content = XMLParser.extractObjectForTags(representation, "content", String.class );
				} catch( NonParsableException e ) {
					content = null;
				}
				if( XMLParser.hasTag(representation, "ext", null, null) ) {
					ext = (String) XMLParser.extractObjectForTags(representation, "ext");
				}
			}
		}
		
		public boolean equals( Object o ) {
			if(this == o) {
				return true;
			}
			if( o instanceof FileRepresentation ) {
				FileRepresentation fr = (FileRepresentation)o;
								
				String line1, line2;
				Character ignore = null;
				if( ext != null ) {
					switch( ext ) {
						case "gff": case "gff3": ignore='#'; break;
						default: ignore = null;
					}
				}
				int i = 0;
				try {
					BufferedReader r1 = new BufferedReader(content == null ? new FileReader(filename) : new StringReader(compressed ? Compression.unzip(content) : content) );
					BufferedReader r2 = new BufferedReader(fr.content == null ? new FileReader(fr.filename) : new StringReader(fr.compressed ? Compression.unzip(fr.content) : fr.content) );
					do {
						while( (line1=r1.readLine()) != null && line1.length()>0 && ((Character)line1.charAt(0)).equals(ignore) );
						while( (line2=r2.readLine()) != null && line2.length()>0 && ((Character)line2.charAt(0)).equals(ignore) );
						i++;
					} while( line1!= null && line2 != null && line1.equals(line2) );	
					r1.close();
					r2.close();
				} catch( IOException e ) {
					e.printStackTrace();
					return false;
				}
				if( line1!= null || line2!= null ) {
					System.err.println("Files differ at line " + i );
					System.err.println(filename + "\t" + line1);
					System.err.println(fr.filename + "\t" + line2);
				}
				return line1==null && line2==null;
			} else {
				return false;
			}
		}

	}
}

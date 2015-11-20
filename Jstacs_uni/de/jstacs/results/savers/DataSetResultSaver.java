package de.jstacs.results.savers;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.annotation.SplitSequenceAnnotationParser;
import de.jstacs.results.DataSetResult;


/**
 * Class for a {@link ResultSaver} working on {@link DataSetResult}.
 * The result's contents are saved using the {@link DataSet#save(java.io.OutputStream, char, de.jstacs.data.sequences.annotation.SequenceAnnotationParser)} method in FastA format.
 * 
 * @author Jan Grau
 *
 */
public class DataSetResultSaver implements ResultSaver<DataSetResult> {

	/**
	 * Registers this {@link ResultSaver} in the {@link ResultSaverLibrary}
	 */
	public static void register(){
		ResultSaverLibrary.register( DataSetResult.class, new DataSetResultSaver() );
	}
	
	private DataSetResultSaver() {
	}
	
	@Override
	public boolean isAtomic() {
		return true;
	}

	@Override
	public String[] getFileExtensions( DataSetResult result ) {
		return new String[]{"fa"};
	}

	@Override
	public boolean writeOutput( DataSetResult result, File path ) {
		try{
			FileOutputStream fos = new FileOutputStream( path );
			if( (result).getParser() == null ){
				(result).getValue().save(fos,'>',new SplitSequenceAnnotationParser( ":", ";" ) );
			}else{
				(result).getValue().save(fos,'>', result.getParser() );
			}
			fos.close();
			return true;
		}catch(IOException e){
			e.printStackTrace();
			return false;
		}
	}

	@Override
	public boolean writeOutput( DataSetResult result, StringBuffer buf ) {
		try{
			ByteArrayOutputStream baos = new ByteArrayOutputStream();
			if( (result).getParser() == null ){
				(result).getValue().save(baos,'>',new SplitSequenceAnnotationParser( ":", ";" ) );
			}else{
				(result).getValue().save(baos,'>', result.getParser() );
			}
			buf.append( baos.toString() );
			return true;
		}catch(IOException e){
			e.printStackTrace();
			return false;
		}
	}
	
	

}

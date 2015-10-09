package de.jstacs.results.savers;

import java.io.File;
import java.util.HashSet;

import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.ResultSetResult;
import de.jstacs.tools.ToolResult;


public class ResultSetResultSaver implements ResultSaver<ResultSetResult> {

	public static void register(){
		ResultSaverLibrary.register( ResultSetResult.class, new ResultSetResultSaver() );
		ResultSaverLibrary.register( ToolResult.class, new ResultSetResultSaver() );
	}
	
	private ResultSetResultSaver() {
	}

	@Override
	public boolean isAtomic() {
		return false;
	}

	@Override
	public String[] getFileExtensions( ResultSetResult result ) {
		return null;
	}

	@Override
	public boolean writeOutput( ResultSetResult result, File dir ) {
		if(!dir.exists()){
			dir.mkdirs();
		}
		
		if(!dir.isDirectory()){
			return false;
		}
		ResultSet set = result.getRawResult()[0];
		HashSet<String> names = new HashSet<String>();
		boolean wroteAll = true;
		for(int i=0;i<set.getNumberOfResults();i++){
			Result res = set.getResultAt( i );
			ResultSaver saver = ResultSaverLibrary.getSaver( res );
			if(saver != null){
				String filename = res.getName().replaceAll( "[\\s\\:\\/]", "_" );
				String temp = filename;
				if(saver.isAtomic()){
					temp = temp + "." + saver.getFileExtensions( res )[0];
				}
				int j=1;
				while(names.contains( temp ) || new File(dir.getAbsolutePath()+File.separator+temp).exists()){//TODO correct?
					temp = filename+"_"+j;
					if(saver.isAtomic()){
						temp = temp + "." + saver.getFileExtensions( res )[0];
					}
					j++;
				}
				filename = temp;
				names.add( filename );
				
				
				wroteAll &= saver.writeOutput( res, new File(dir.getAbsolutePath()+File.separator+filename) );
				
			}
		}
		return wroteAll;
	}

	@Override
	public boolean writeOutput( ResultSetResult result, StringBuffer buf ) {
		throw new RuntimeException( "Not possible" );
	}

	
	
}

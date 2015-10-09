package de.jstacs.results.savers;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

import de.jstacs.results.Result;


public class ResultSaverLibrary {

	static{
		TextResultSaver.register();
		ResultSetResultSaver.register();
		ListResultSaver.register();
		DataSetResultSaver.register();
		PlotGeneratorResultSaver.register();
	}
	
	private static HashMap<Class<? extends Result>, ResultSaver> map;
	
	
	public static <T extends Result> void register(Class<? extends T> clazz, ResultSaver<T> renderer){
		if(map == null){
			map = new HashMap<Class<? extends Result>, ResultSaver>();
		}
		map.put( clazz, renderer );
	}
	
	public static <T extends Result> ResultSaver<T> getSaver(T result){
		
		ResultSaver ren = map.get( result.getClass() );
		if(ren == null){
			Set<Class<? extends Result>> clazzes = map.keySet();
			Iterator<Class<? extends Result>> it = clazzes.iterator();
			while(it.hasNext()){
				Class<? extends Result> clazz = it.next();
				if(clazz.isAssignableFrom( result.getClass() )){
					return map.get( clazz );
				}
			}
		}
		if(ren == null){
			System.err.println( "Did not find a saver for "+result.getClass()+". Custom savers need to be registered by "+ResultSaverLibrary.class.getName()+".register()." );
		}
		return ren;
	}

}

package de.jstacs.results.savers;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

import de.jstacs.results.Result;

/**
 * Library of all registered {@link ResultSaver}s.
 * 
 * @author Jan Grau
 *
 */
public class ResultSaverLibrary {

	//registers standard ResultSavers
	static{
		TextResultSaver.register();
		ResultSetResultSaver.register();
		ListResultSaver.register();
		DataSetResultSaver.register();
		PlotGeneratorResultSaver.register();
	}
	
	private static HashMap<Class<? extends Result>, ResultSaver> map;
	
	/**
	 * Registers the supplied {@link ResultSaver} for the given class in the library.
	 * @param clazz the class
	 * @param saver the {@link ResultSaver}
	 */
	public static <T extends Result> void register(Class<? extends T> clazz, ResultSaver<T> saver){
		if(map == null){
			map = new HashMap<Class<? extends Result>, ResultSaver>();
		}
		map.put( clazz, saver );
	}
	
	/**
	 * Returns the most suitable {@link ResultSaver} (if any) currently registered in the library.
	 * If a {@link ResultSaver} for the given {@link Result} type is registered, this one is returned.
	 * Otherwise, the list of {@link ResultSaver}s is searched for a {@link ResultSaver} registered for
	 * a superclass of <code>T</code> and the first hit is returned.
	 * 
	 * @param result the result to be saved
	 * @return the appropriate {@link ResultSaver} or <code>null</code> if not such {@link ResultSaver} has been registered
	 */
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

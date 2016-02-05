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

package de.jstacs.utils;

import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.lang.reflect.Constructor;
import java.lang.reflect.Modifier;
import java.net.JarURLConnection;
import java.net.MalformedURLException;
import java.net.URI;
import java.net.URISyntaxException;
import java.net.URL;
import java.net.URLEncoder;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.jar.JarEntry;
import java.util.jar.JarFile;

import de.jstacs.InstantiableFromParameterSet;
import de.jstacs.Singleton;
import de.jstacs.Singleton.SingletonHandler;
import de.jstacs.io.RegExFilenameFilter;
import de.jstacs.io.RegExFilenameFilter.Directory;
import de.jstacs.parameters.InstanceParameterSet;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SelectionParameter;

/**
 * Utility-class with static methods to
 * <ul>
 * <li>find all sub-classes of a certain class (or interface) within the scope
 * of the current class-loader</li>
 * <li>find all sub-classes of a certain class (or interface) within the scope
 * of the current class-loader that can be instantiated, i.e. that are neither
 * interfaces nor abstract</li>
 * <li>filter a set of classes by inheritance from a super-class</li>
 * <li>obtain the class of an {@link InstanceParameterSet} that can be used to
 * instantiate a sub-class of {@link InstantiableFromParameterSet}.</li>
 * <li>obtain a {@link SelectionParameter} using all possible
 * {@link ParameterSet}s (for classes that are a subclass of a specified
 * superclass) as elements</li>
 * </ul>
 * 
 * The methods may fail to find a certain sub-class, if
 * <ul>
 * <li>it was loaded by another {@link ClassLoader} than the caller</li>
 * <li>it was loaded from a physically other place than the super-class, e.g.
 * another jar-file.
 * </ul>
 * 
 * @author Jens Keilwagen, Jan Grau
 */
public class SubclassFinder {

	/**
	 * This field can be set to include a path into the search performed in {@link #findSubclasses(Class, String)}
	 * thereby enabling to find self-implemented classes not included in the Jstacs class hierarchy.  
	 * The default value of this field is <code>null</code>.
	 * 
	 * <br><br>
	 * 
	 * <b>
	 * Please use this field very careful, i.e., after setting a specific value and calling 
	 * {@link #findSubclasses(Class, String)} it is highly recommended to reset the value. If this is not done, preceeding searches might take longer.
	 * </b>
	 */
	public static String includePath = null;
	
	/**
	 * Returns all sub-classes of <code>T</code> that can be instantiated, i.e.
	 * are neither an interface nor abstract, and that are located in a package
	 * below <code>startPackage</code>.
	 * 
	 * @param <T>
	 *            The class to obtain the sub-classes for
	 * @param clazz
	 *            the {@link Class} object for T
	 * @param startPackage
	 *            the package under which to search
	 * @return the {@link Class} objects for the sub-classes
	 * @throws ClassNotFoundException
	 *             if one of the classes is present in the file system or jar
	 *             but cannot be loaded by the class loader
	 * @throws IOException
	 *             is thrown if the classes are searched for in a jar file, but
	 *             that file could not be accessed or read
	 */
	public static <T> LinkedList<Class<? extends T>> findInstantiableSubclasses( Class<T> clazz, String startPackage ) throws ClassNotFoundException,
			IOException {
		LinkedList<Class<? extends T>> list = findSubclasses( clazz, startPackage );
		LinkedList<Class<? extends T>> list2 = new LinkedList<Class<? extends T>>();
		Iterator<Class<? extends T>> it = list.iterator();
		while( it.hasNext() ) {
			Class<? extends T> c = it.next();
			if( !c.isInterface() && !Modifier.isAbstract( c.getModifiers() ) ) {
				list2.add( c );
			}
		}
		return list2;
	}

	/**
	 * Returns a {@link LinkedList} of the {@link Class} objects for all classes
	 * in <code>listToFilter</code> that are sub-classes of
	 * <code>superClass</code>.
	 * 
	 * @param <S>
	 *            the class that is used as filter
	 * @param <T>
	 *            a common super-class of all classes in
	 *            <code>listToFilter</code>
	 * @param superclass
	 *            the additional class to use as a filter criterion
	 * @param listToFilter
	 *            the list of classes
	 * @return the filtered list
	 */
	public static <S, T> LinkedList<Class<? extends S>> filterBySuperclass( Class<S> superclass, LinkedList<Class<? extends T>> listToFilter ) {

		LinkedList<Class<? extends S>> list2 = new LinkedList<Class<? extends S>>();
		Iterator<Class<? extends T>> it = listToFilter.iterator();
		while( it.hasNext() ) {
			Class<? extends T> c = it.next();
			if( superclass.isAssignableFrom( c ) ) {
				list2.add( (Class<? extends S>)c );
			}
		}
		return list2;
	}

	/**
	 * Returns a {@link LinkedList} of the classes of the
	 * {@link InstanceParameterSet}s that can be used to instantiate the
	 * sub-class of {@link InstantiableFromParameterSet} that is given by
	 * <code>clazz</code>
	 * 
	 * @param clazz
	 *            the {@link Class} object of the sub-class of
	 *            {@link InstantiableFromParameterSet}
	 * @return a {@link LinkedList} of {@link Class}es of the corresponding
	 *         {@link InstanceParameterSet}s
	 */
	public static LinkedList<Class<? extends InstanceParameterSet>> getParameterSetFor( Class<? extends InstantiableFromParameterSet> clazz ) {
		Constructor[] cons = clazz.getConstructors();
		Class[] types = null;
		LinkedList<Class<? extends InstanceParameterSet>> list = new LinkedList<Class<? extends InstanceParameterSet>>();
		boolean add = false;
		for( int i = 0; i < cons.length; i++ ) {
			if( ( types = cons[i].getParameterTypes() ).length == 1 ) {
				if( InstanceParameterSet.class.isAssignableFrom( types[0] ) ) {
					list.add( types[0] );
					add = true;
				}
			}
		}
		if( !add && Singleton.class.isAssignableFrom(clazz) ) {
			try {
				InstantiableFromParameterSet ifps = (InstantiableFromParameterSet)SingletonHandler.getSingelton((Class<? extends Singleton>) clazz);
				list.add( ifps.getCurrentParameterSet().getClass() );
			} catch ( Exception e ) {
				throw new RuntimeException( e );
			}
		}
		return list;
	}

	/**
	 * Returns all sub-classes of <code>T</code> including interfaces and
	 * abstract classes that are located in a package below
	 * <code>startPackage</code>.
	 * 
	 * @param <T>
	 *            The class to obtain the sub-classes for
	 * @param clazz
	 *            the {@link Class} object for T
	 * @param startPackage
	 *            the package under which to search
	 * @return the {@link Class} objects for the sub-classes
	 * @throws ClassNotFoundException
	 *             if one of the classes is present in the file system or jar
	 *             but cannot be loaded by the class loader
	 * @throws IOException
	 *             is thrown if the classes are searched for in a jar file, but
	 *             that file could not be accessed or read 
	 */
	public static <T> LinkedList<Class<? extends T>> findSubclasses( Class<T> clazz, String startPackage ) throws ClassNotFoundException, IOException {
		HashSet<Class<? extends T>> hash = new HashSet<Class<? extends T>>();
		
		find( clazz, startPackage, hash );
		find( clazz, includePath, hash );

		return new LinkedList<Class<? extends T>>( hash );
	}
	
	private static <T> void find( Class<T> clazz, String startPackage, HashSet<Class<? extends T>> hash ) throws ClassNotFoundException, IOException {
		if( startPackage == null )  {
			return;
		}
		String name = startPackage;
		if( !name.startsWith( "/" ) ) {
			name = "/" + name;
		}
		name = name.replace( ".", "/" );
		URL[] urls = getResources( name );
		for(URL url : urls){
			
			if( url != null ) {
				String s = url.toString();
				if( s.startsWith("file:") ) {
					s = s.substring(5);
				}
				File dir = new File( s );
				
				if( dir.exists() ) {
					File[] files = dir.listFiles();
					for( int i = 0; i < files.length; i++ ) {
						if( files[i].isDirectory() ) {
							//System.out.println(startPackage+"."+files[i].getName());
							find( clazz, startPackage + "." + files[i].getName(), hash );
						} else if( files[i].isFile() && files[i].getName().endsWith( ".class" ) ) {
							add( clazz, hash, startPackage + "." + files[i].getName().substring( 0, files[i].getName().lastIndexOf( "." ) ) );
						}
					}
				} else {
					JarURLConnection con = (JarURLConnection)url.openConnection();
					JarFile jar = con.getJarFile();
					String starts = con.getEntryName();
					Enumeration<JarEntry> en = jar.entries();
					while( en.hasMoreElements() ) {
						JarEntry entry = en.nextElement();
						String entryname = entry.getName();
						if( (starts == null || entryname.startsWith( starts ) ) && entryname.endsWith( ".class" ) ) {
							String classname = entryname.substring( 0, entryname.length() - 6 );
							if( classname.startsWith( "/" ) ) {
								classname = classname.substring( 1 );
							}
							add( clazz, hash, classname.replace( "/", "." ) );
						}
						
					}
				}
			}
		}
	}

	@SuppressWarnings( "unchecked" )
	private static <T> void add( Class<T> clazz, HashSet<Class<? extends T>> list, String className ) {
		try {
			Class c = Class.forName( className );
			if( clazz.isAssignableFrom( c ) ) {
				list.add( c );
			}
		} catch( Throwable t ) {
			/*
			System.out.println( clazz + "\t" + className );
			t.printStackTrace();
			throw new RuntimeException( t.getMessage() );
			*/
		}
	}
	
	private static URL[] getResources(String startPackage) throws MalformedURLException{
		String cp = System.getProperty( "java.class.path" );
		String ext = System.getProperty( "java.ext.dirs" );
		String sun = System.getProperty( "sun.boot.class.path" );
		
		LinkedList<URL> list = new LinkedList<URL>();
		
		String startStart = startPackage;
		if(startStart.startsWith( "/" )){
			startStart = startStart.substring( 1 );
		}
		
		if(sun != null){
			cp = cp+System.getProperty( "path.separator" )+sun;
		}
		
		URL url;
		JarURLConnection con;
		
		if(ext != null){
			String[] extDirs = ext.split( System.getProperty( "path.separator" ) );
			for(String dir : extDirs){
				File temp = new File(dir);
				if(temp.isDirectory()){
					File[] jars = temp.listFiles( (FileFilter) new RegExFilenameFilter( "", Directory.FORBIDDEN, true, ".*(\\.jar|\\.zip)" ) );
					for(File jar : jars){
						url = new URL("jar:file:"+jar.getAbsolutePath()+"!"+startPackage);
						try{
							con = (JarURLConnection)url.openConnection();
							con.getJarFile();
							list.add( url );
						}catch(IOException e){
							
						}
					}
				}
			}
		}
		
		
		
		String[] cps = cp.split( System.getProperty( "path.separator" ) );
		File curr;
		outerloop:
		for(int i=0;i<cps.length;i++){
			cps[i] = cps[i].trim();
			curr = new File( cps[i] );
			if( curr.exists() && curr.isDirectory() && (new File(cps[i]+System.getProperty( "file.separator" )+startPackage).exists()) ){
				list.add( new URL( "file:"+cps[i]+startPackage ) );
			}else if(cps[i].endsWith( ".jar" ) || cps[i].endsWith( ".zip" )){
				
				url = new URL( "jar:file:"+cps[i]+"!"+startPackage );
				try{
					con = (JarURLConnection)url.openConnection();
					con.getJarFile();
					list.add( url );
				}catch(IOException e){
					try{
						url = new URL( "jar:file:"+cps[i]+"!/" );
						con = (JarURLConnection)url.openConnection();
						JarFile fj = con.getJarFile();
						Enumeration<JarEntry> en = fj.entries();
						while(en.hasMoreElements()){
							if(en.nextElement().getName().startsWith( startStart )){
								list.add( url );
								continue outerloop;
							}
						}
					}catch(Exception ex){
					}
					
					
				}
			}
		}
		
		return list.toArray( new URL[0] );
		
	}

	/**
	 * This method creates an {@link SelectionParameter} that contains
	 * {@link de.jstacs.parameters.InstanceParameterSet} for each possible
	 * class. The classes are specified by
	 * {@link SubclassFinder#findInstantiableSubclasses(Class, String)} and
	 * {@link SubclassFinder#filterBySuperclass(Class, LinkedList)}.
	 * 
	 * @param <T>
	 *            The class to use the sub-classes in the
	 *            {@link SelectionParameter}
	 * @param clazz
	 *            the {@link Class} object for <code>T</code>
	 * @param startPackage
	 *            the package under which to start the search
	 * @param name
	 *            the name of the {@link SelectionParameter}
	 * @param comment
	 *            the comment for the {@link SelectionParameter}
	 * @param required
	 *            whether the {@link SelectionParameter} is required
	 * 
	 * @return a {@link SelectionParameter} that contains
	 *         {@link de.jstacs.parameters.InstanceParameterSet} for each
	 *         possible class.
	 *         
	 * @throws Exception
	 *             if the {@link SelectionParameter} could not be created properly
	 * 
	 * @see SubclassFinder#findInstantiableSubclasses(Class, String)
	 * @see SubclassFinder#filterBySuperclass(Class, LinkedList)
	 * @see #getInstanceParameterSets(Class, String)
	 */
	@SuppressWarnings("unchecked")
	public static <T> SelectionParameter getSelectionParameter( Class<? extends ParameterSet> clazz, String startPackage, String name, String comment, boolean required ) throws Exception {
		LinkedList<?> list = SubclassFinder.findInstantiableSubclasses( clazz, startPackage );
		Class<? extends ParameterSet>[] classes = new Class[list.size()];
		Iterator<?> it = list.iterator();
		for( int i = 0; i < classes.length; i++ ){
			classes[i] = (Class<? extends ParameterSet>)it.next();
		}
		Arrays.sort( classes, ClassNameComparator.DEFAULT );
		return new SelectionParameter( name, comment, required, classes );
	}
	
	private static class ClassNameComparator implements Comparator<Class> {

		private static ClassNameComparator DEFAULT = new ClassNameComparator();
		
		private ClassNameComparator(){}
		
		@Override
		public int compare(Class o1, Class o2) {
			return o1.getSimpleName().compareTo( o2.getSimpleName() );
		}
	}
	
	/**
	 * This method returns a list of {@link InstanceParameterSet}s that can be used to create a subclass of <code>clazz</code>.
	 * The classes are specified by
	 * {@link SubclassFinder#findInstantiableSubclasses(Class, String)} and
	 * {@link SubclassFinder#filterBySuperclass(Class, LinkedList)}.
	 * 
	 * @param <T>
	 *            The class to use the sub-classes in the
	 *            {@link SelectionParameter}
	 * @param clazz
	 *            the {@link Class} object for <code>T</code>
	 * @param startPackage
	 *            the package under which to start the search
	 * 
	 * @return a {@link LinkedList} that contains
	 *         {@link de.jstacs.parameters.InstanceParameterSet} for each
	 *         possible class.
	 * 
	 * @throws InstantiationException
	 *             if any {@link InstanceParameterSet} has no nullary
	 *             constructor; or if the instantiation fails for some other
	 *             reason
	 * @throws IllegalAccessException
	 *             if any {@link InstanceParameterSet} or its nullary
	 *             constructor is not accessible
	 * @throws ClassNotFoundException
	 *             if one of the classes is present in the file system or jar
	 *             but cannot be loaded by the class loader
	 * @throws IOException
	 *             if the classes are searched for in a jar file, but that file
	 *             could not be accessed or read
	 * 
	 * @see SubclassFinder#findInstantiableSubclasses(Class, String)
	 * @see SubclassFinder#filterBySuperclass(Class, LinkedList)
	 */
	public static <T> LinkedList<InstanceParameterSet<? extends T>> getInstanceParameterSets( Class<T> clazz, String startPackage ) throws ClassNotFoundException, IOException, InstantiationException, IllegalAccessException {
		LinkedList<Class<? extends T>> classes = SubclassFinder.findInstantiableSubclasses( clazz, startPackage );
		LinkedList<Class<? extends InstantiableFromParameterSet>> filteredClasses = SubclassFinder.filterBySuperclass( InstantiableFromParameterSet.class,
				classes );
		LinkedList<InstanceParameterSet<? extends T>> sets = new LinkedList<InstanceParameterSet<? extends T>>();
		Iterator<Class<? extends InstantiableFromParameterSet>> it = filteredClasses.iterator();
		Iterator<Class<? extends InstanceParameterSet>> psIt;
		Class<? extends InstanceParameterSet> c;
		while( it.hasNext() ) {
			psIt = SubclassFinder.getParameterSetFor( it.next() ).iterator();
			while( psIt.hasNext() ) {
				c = psIt.next();
				if( Singleton.class.isAssignableFrom(c) ) {
					try {
						sets.add( (InstanceParameterSet) SingletonHandler.getSingelton( (Class<? extends Singleton>) c ) );
					} catch ( Exception e ) {
						throw new RuntimeException( e );
					}
				} else {
					sets.add( c.newInstance() );
				}
			}
		}
		return sets;
	}
}

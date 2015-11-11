package de.jstacs.data;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;

import de.jstacs.data.OboParser.GONode;


public class OboParser {

	private HashMap<String, GONode> map;
	
	public OboParser(String filename) throws IOException{
		
		BufferedReader read = new BufferedReader( new FileReader( filename ) );
		
		String str = null;
		
		String go = null;
		String name = null;
		
		map = new HashMap<String, OboParser.GONode>();
		
		LinkedList<GONode> isas = new LinkedList<OboParser.GONode>();
		
		while( (str = read.readLine()) != null ){
			str = str.trim();
			if(str.length() == 0){
				if(go != null && go.startsWith( "GO:" )){
					
					if(map.containsKey( go )){
						map.get( go ).setName( name );
						map.get( go ).setParents( isas.toArray( new GONode[0] ) );
					}else{
						map.put( go, new GONode( go, name, isas.toArray( new GONode[0] ) ) );
					}
					
				}
				
				
				go = null;
				name = null;
				isas.clear();
			}else if(str.startsWith( "id:" )){
				go = str.substring( 4 );
			}else if(str.startsWith( "name:" )){
				name = str.substring( 6 );
			}else if(str.startsWith( "is_a:" )){
				String temp = str.substring( 6 );
				temp = temp.substring( 0,temp.indexOf( ' ' ));
				
				if(map.containsKey( temp )){
					isas.add( map.get( temp ) );
				}else{
					GONode node = new GONode( temp );
					isas.add( node );
					map.put( temp, node );
				}
				
			}else if(str.startsWith( "is_obsolete: true" )){
				go = null;
				name = null;
				isas.clear();
			}
		}
		
		
		read.close();
		
		
		Iterator<String> it = map.keySet().iterator();
		
		int i=1;
		
		while(it.hasNext()){
			GONode node = map.get( it.next() );
			
			if(node.isRoot()){
				System.out.println(node.go);
				node.propagateIds(i+"",i+"");
				i++;
			}
			
			
		}
		
		
	}
	
	
	
	public static class GONode{
		
		private final String go;
		private String name;
		private HashSet<String> ids;
		
		private GONode[] parents;
		private HashSet<GONode> children;
		private LinkedList<String> lastIdPart;

		public GONode(String go){
			this.go = go;
			this.children = new HashSet<OboParser.GONode>();
			this.ids = new HashSet<String>();
			this.lastIdPart = new LinkedList<String>();
		}
		

		public GONode(String go, String name, GONode[] parents){
			this(go);
			this.name = name;
			setParents( parents );
		}
		
		public String toString(){
			return go+"\n"+name+"\n"+Arrays.toString( ids.toArray( new String[0] ) );
		}
		
		public String getGO(){
			return go;
		}
		
		public String[] getHierarchyNames(){
			if(parents == null || parents.length == 0){
				if(name == null){
					System.out.println("null: "+go);
				}
				return new String[]{name};
			}
			
			LinkedList<String> names = new LinkedList<String>();
			
			for(int i=0;i<parents.length;i++){
				String[] par = parents[i].getHierarchyNames();
				for(int j=0;j<par.length;j++){
					names.add( par[j]+"."+name );
				}
			}
			return names.toArray( new String[0] );
		}
		
		public void propagateIds( String id, String last ) {
			this.lastIdPart.add(last);
			this.ids.add( id );
			Iterator<GONode> it = children.iterator();
			int i=1;
			while(it.hasNext()){
				GONode curr = it.next();
				curr.propagateIds( id+"."+i, i+"" );
				i++;
			}
		}
		
		public LinkedList<String> getLastIDPart(){
			return lastIdPart;
		}
		
		public String[] getIDs(){
			return ids.toArray( new String[0] );
		}
		
		public String[] getMajorIDAndName(){
			String[][] ids = new String[this.ids.size()][];
			
			Iterator<String> it = this.ids.iterator();
			int i=0;
			while(it.hasNext()){
				String curr = it.next();
				String[] parts = curr.split( "\\." );
				ids[i] = parts;
				i++;
			}
			
			String[] dom = findDominant(ids,0);
			String ret = dom[0];
			for(i=1;i<dom.length;i++){
				ret += "."+dom[i];
			}
			
			int lv = dom.length-1;
			
			String name = this.getNameFor( dom, lv );
			
			return new String[]{ret,name};
		}
		
		public String getNameFor(String[] dom, int lv){
			
			String name = this.name;
			if(lv<=0 || parents == null || parents.length==0){
				return name;
			}else{
				String partDom = dom[0];
				for(int i=1;i<lv;i++){
					partDom += "."+dom[i];
				}
				for(int i=0;i<parents.length;i++){
					if(parents[i].ids.contains( partDom ) ){
						return parents[i].getNameFor( dom, lv-1 )+"."+name;
					}
				}
				throw new RuntimeException(dom[lv-1]+" "+lv+" "+Arrays.toString( dom )+" <-> "+this.ids+" "+partDom);
			}
		}
		
		private static String[] findDominant( String[][] ids2, int level ) {
			HashMap<String, Integer> counters = new HashMap<String, Integer>();
			int max = 0;
			for(int i=0;i<ids2.length;i++){
				if(level < ids2[i].length){
					if(counters.containsKey( ids2[i][level] )){
						counters.put( ids2[i][level], counters.get( ids2[i][level] )+1 );
					}else{
						counters.put( ids2[i][level],1);
					}
					if(counters.get( ids2[i][level] ) > max){
						max = counters.get( ids2[i][level] );
					}
				}
			}
			if(counters.size() == 0){
				return new String[0];
			}
			
			Iterator<String> it = counters.keySet().iterator();
			
			String best = null;
			
			while(it.hasNext()){
				String key = it.next();
				if( counters.get( key ) == max ){
					best = key;
					break;
				}
			}
			
			LinkedList<String[]> lili = new LinkedList<String[]>();
			
			for(int i=0;i<ids2.length;i++){
				if(level < ids2[i].length){
					if(best.equals( ids2[i][level] )){
						lili.add( ids2[i] );
					}
				}
			}
			
			String[] temp = findDominant( lili.toArray( new String[0][0] ), level+1 );
			String[] next = new String[temp.length+1];
			next[0] = best;
			System.arraycopy( temp, 0, next, 1, temp.length );
			return next;
		}


		public boolean isRoot(){
			return parents == null || parents.length == 0;
		}
		
		
		public String getName() {
			return name;
		}

		
		public void setName( String name ) {
			this.name = name;
		}

		
		public GONode[] getParents() {
			return parents;
		}

		
		public void setParents( GONode[] parents ) {
			this.parents = parents;
			if(parents != null){
				for(int i=0;i<this.parents.length;i++){
					this.parents[i].addChild(this);
				}
			}
		}

		private void addChild( GONode goNode ) {
			this.children.add( goNode );
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + ( ( go == null ) ? 0 : go.hashCode() );
			return result;
		}

		@Override
		public boolean equals( Object obj ) {
			if( this == obj ) return true;
			if( obj == null ) return false;
			if( getClass() != obj.getClass() ) return false;
			GONode other = (GONode)obj;
			if( go == null ) {
				if( other.go != null ) return false;
			} else if( !go.equals( other.go ) ) return false;
			return true;
		}
		
		
		
		
	}



	public GONode getNodeFor( String string ) {
		return map.get( string );
	}
	
	
}

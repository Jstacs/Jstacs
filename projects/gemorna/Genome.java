package projects.gemorna;

import java.util.HashMap;

import projects.gemoma.Tools;

public class Genome {

	public static Genome genome;
	public static HashMap<String,Character> code;
	public static boolean[][][] isStart;
	public static boolean[][][] isStop;
	
	private HashMap<String,char[]> chromosomes;
	
	public static void init(String path) throws Exception {
		genome = new Genome(path);
		code = Tools.getCode(Genome.class.getClassLoader().getResourceAsStream( "projects/gemorna/genetic_code.txt"));
		code.put("NNN", 'X');
		

		int maxIdx = 256;
		
		isStart = new boolean[maxIdx][maxIdx][maxIdx];
		isStop = new boolean[maxIdx][maxIdx][maxIdx];
		
		
		

		
		for(String key : code.keySet()) {
			Character val = code.get(key);
			if(val == '*') {
				char[] keys = key.toCharArray();
				isStop[keys[0]][keys[1]][keys[2]] = true;
			}else if(val == 'M') {
				char[] keys = key.toCharArray();
				isStart[keys[0]][keys[1]][keys[2]] = true;
			}
		}
		
	}
	
	public String[] getChromosomeNames() {
		return this.chromosomes.keySet().toArray(new String[0]);
	}
	
	private Genome(String path) throws Exception {
		HashMap<String,String> temp = Tools.getFasta(path, 5, ".*");
		this.chromosomes = new HashMap<String,char[]>();
		for(String key : temp.keySet()) {
			this.chromosomes.put(key, temp.get(key).toCharArray());
		}
	}
	
	
	public char[] getChromosome(String id) {
		return chromosomes.get(id);
	}
	
}

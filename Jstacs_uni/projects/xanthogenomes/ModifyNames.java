package projects.xanthogenomes;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

import de.jstacs.io.FileManager;
import de.jstacs.io.NonParsableException;
import de.jstacs.parameters.FileParameter;

public class ModifyNames {

	public static void main(String[] args) throws Exception {
		TALEFamilyBuilder builder = new TALEFamilyBuilder(FileManager.readFile(args[0]));
		
		BufferedReader read = new BufferedReader(new FileReader(args[1]));
		
		HashMap<String,String[]> map = new HashMap<>();
		String str = null;
		while( (str = read.readLine()) != null ){
			String[] parts = str.split("\t");
			map.put(parts[0], new String[]{parts[1],parts[2]});
		}
		
		
		TALE[] tales = builder.getAllTALEs();
		
		for(int i=0;i<tales.length;i++){
			String id = tales[i].getId();
			String acc = tales[i].getAccession();
			String strain = tales[i].getStrain();
			if(map.keySet().contains(strain)){
				id = id.replaceAll(" "+strain+" ", " "+map.get(strain)[0]+" ");
				acc = map.get(strain)[1];
				strain = map.get(strain)[0];
	
				//System.out.println(id+"###"+strain+"###"+acc);
				tales[i].setAccession(acc);
				tales[i].setId(id);
				tales[i].setStrain(strain);
			}
		}

		FileManager.writeFile(args[0]+"_mod.xml", builder.toXML());
		
	}

}

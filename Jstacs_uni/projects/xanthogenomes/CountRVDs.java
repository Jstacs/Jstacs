package projects.xanthogenomes;

import de.jstacs.io.FileManager;
import projects.xanthogenomes.TALE.Type;

public class CountRVDs {

	public static void main(String[] args) throws Exception {
		
		TALEFamilyBuilder builder = new TALEFamilyBuilder(FileManager.readFile(args[0]));
		
		TALE[] tales = builder.getAllTALEs();
		
		for(int i=0;i<tales.length;i++){
			int num = 0;
			int numLong = 0;
			int numShort = -1;
			for(int j=0;j<tales[i].getNumberOfRepeats();j++){
				if(tales[i].getRepeat(j).getType() == Type.LONG){
					numLong++;
				}else if(tales[i].getRepeat(j).getType() == Type.SHORT){
					numShort++;
				}
				num++;
			}
			System.out.println(tales[i].getId()+"\t"+num+"\t"+numLong+"\t"+numShort);
		}

	}

}

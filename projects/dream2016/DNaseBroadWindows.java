package projects.dream2016;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import de.jstacs.utils.ToolBox;

public class DNaseBroadWindows {

	public static void main(String[] args) throws IOException {
		
		GZIPInputStream stream = new GZIPInputStream(new FileInputStream(args[0]));
		BufferedReader reader = new BufferedReader(new InputStreamReader(stream));

		GZIPOutputStream os = new GZIPOutputStream(new FileOutputStream(args[0]+"_bw"+args[1]+".gz"));
		PrintWriter wr = new PrintWriter(os);
		
		int bins = Integer.parseInt(args[1]);
		
		double[] minProfile = new double[2*bins+1];
		double[] maxProfile = new double[2*bins+1];
		
		int end = 0;
		
		String str = null;
		while( true ){
			str = reader.readLine();
			if(str == null || str.startsWith("[")){
				if(end>0){
					end--;
				}
				for(;end>bins;end--){
					double minVal = ToolBox.min(minProfile.length-end, minProfile.length, minProfile);
					double maxVal = ToolBox.max(maxProfile.length-end, maxProfile.length, maxProfile);
					
					double minBefore = ToolBox.min(minProfile.length-end, minProfile.length-end+bins, minProfile);
					double minAfter = 0;
					if(0 > -end+bins+1){
						minAfter = ToolBox.min( minProfile.length-end+bins+1, minProfile.length, minProfile);
					}
					
					double maxBefore = ToolBox.max(maxProfile.length-end, maxProfile.length-end+bins, maxProfile);
					double maxAfter = 0;
					if(0 > -end+bins+1){
						maxAfter = ToolBox.max( maxProfile.length-end+bins+1, maxProfile.length, maxProfile);
					}
					
					wr.println(minVal+"\t"+maxVal+"\t"+minBefore+"\t"+minAfter+"\t"+maxBefore+"\t"+maxAfter);
				}
				if(str == null){
					break;
				}else{
					wr.println(str);
					Arrays.fill(minProfile, 0);
					Arrays.fill(maxProfile, 0);
					end = 0;
				}
			}else{
				String[] parts = str.split("\t");
				double min = Double.parseDouble(parts[0]);
				double max = Double.parseDouble(parts[2]);
				if(end < minProfile.length){
					minProfile[end] = min;
					maxProfile[end] = max;
					end++;
					
				}else{
					for(int i=0;i<minProfile.length-1;i++){
						minProfile[i] = minProfile[i+1];
						maxProfile[i] = maxProfile[i+1];
					}
					minProfile[minProfile.length-1] = min;
					maxProfile[maxProfile.length-1] = max;
				}
				if(end>bins){
					double minVal = ToolBox.min(0,end,minProfile);
					double maxVal = ToolBox.max(0,end,maxProfile);
					double minBefore = 0;
					double minAfter = ToolBox.min(Math.min(bins,end-bins-1)+1,end,minProfile);;
					if(end > bins+1){
						minBefore = ToolBox.min(0,Math.min(bins,end-bins-1),minProfile);
					}
					double maxBefore = 0;
					double maxAfter = ToolBox.max(Math.min(bins,end-bins-1)+1,end,maxProfile);;
					if(end > bins+1){
						maxBefore = ToolBox.max(0,Math.min(bins,end-bins-1),maxProfile);
					}
					wr.println(minVal+"\t"+maxVal+"\t"+minBefore+"\t"+minAfter+"\t"+maxBefore+"\t"+maxAfter);
				}
			}
		}
		
		reader.close();
		wr.close();
	}

}

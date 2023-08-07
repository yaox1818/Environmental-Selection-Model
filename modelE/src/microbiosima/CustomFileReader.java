/*******************************************************************************
 *
 * Copyright (C) 2014, 2015 Qinglong Zeng, Jeet Sukumaran, Steven Wu and Allen Rodrigo
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package microbiosima;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/*
 * @author John
 */
public class CustomFileReader {
	private String[] commands;
        public int numberOfLines;

	public CustomFileReader(String file_path, int number_of_line) throws IOException {
		BufferedReader bf = new BufferedReader(new FileReader(file_path));
		commands = new String[number_of_line];
                numberOfLines=number_of_line;
		int i = 0;
		while (i<number_of_line) {
			commands[i] = bf.readLine();
			i++;
		}
		bf.close();
	}
        
        public CustomFileReader(String file_path) throws IOException {
		BufferedReader bf = new BufferedReader(new FileReader(file_path));
		String aLine;
                List<String> commandList=new ArrayList<>();
		int i = 0;
		while ((aLine = bf.readLine()) != null) {
			commandList.add(aLine);
			i++;
		}
                commands = new String[i];
                commandList.toArray(commands);
                numberOfLines=i;
		bf.close();
	}
        
	/**
	 * @return the commands
	 */
	public String getCommand(int index) {
		return commands[index];
	}
	
	public String[] getCommandSplit(int index){
		return commands[index].split("\t");
	}
        
        public double[] getNumericArray(int index){
            String[] stringArray=getCommandSplit(index);
            double[] numericArray=new double[stringArray.length];
            for (int i=0;i<stringArray.length;i++){
                numericArray[i]=Double.parseDouble(stringArray[i]);
            }
            return numericArray;
        }
	
	

}

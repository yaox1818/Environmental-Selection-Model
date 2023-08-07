package microbiosima;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.Arrays;

import org.apache.commons.cli.*;
import utils.random.MathUtil;

import org.apache.commons.cli.Option.Builder;


public class TestE extends Microbiosima {

    private static final String VERSION = "2.0";

    /**
     * @param args
     *            the command line arguments
     * @throws java.io.FileNotFoundException
     * @throws java.io.UnsupportedEncodingException
     */
//修改过k的取值范围

    /**
     将k的过程修改成了一个方法，并且不用每次进行计算，直接继承的k之后的fitness
     */
    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
        int populationSize=5000;//Integer.parseInt(parameters[1]);5000
        int microSize=1000000000;//Integer.parseInt(parameters[2]);1000000000
        int numberOfSpecies=150;//Integer.parseInt(parameters[3]);
        int numberOfGeneration=20000;//200000
        int Ngene=25;//总共有25种可用trait
        int numberOfObservation=200;//5000
        int numberOfReplication=5;//5
        double Ngenepm=5;//每种OTU微生物有5种trait
        double pctEnv=0;//x=1-pctEnv,环境贡献
        double pctPool=0;//y 亲代对环境的贡献
        double msCoeffInHost=1;
        double msCoeffInEnv=1;
        double hsCoeff=1;
        boolean HMS_or_TMS=false;



        Options options = new Options();

        Option help = new Option("h", "help", false, "print this message");
        Option version = new Option("v", "version", false,
                "print the version information and exit");
        options.addOption(help);
        options.addOption(version);

        options.addOption(Option.builder("o").longOpt("obs").hasArg()
                .argName("OBS").desc("Number generation for observation [default: 100]")
                .build());
        options.addOption(Option.builder("r").longOpt("rep").hasArg()
                .argName("REP").desc("Number of replication [default: 1]")
                .build());

        Builder C = Option.builder("c").longOpt("config")
                .numberOfArgs(6).argName("Pop Micro Spec Gen")
                .desc("Four Parameters in the following orders: "
                        + "(1) population size, (2) microbe size, (3) number of species, (4) number of generation, (5) number of total traits, (6)number of traits per microbe"
                        + " [default: 500 1000 150 10000 10 5]");
        options.addOption(C.build());

        HelpFormatter formatter = new HelpFormatter();
        String syntax = "microbiosima pctEnv pctPool";
        String header = "\nSimulates the evolutionary and ecological dynamics of microbiomes within a population of hosts.\n\n"+
                "required arguments:\n"+"  pctEnv             Percentage of environmental acquisition\n"+
                "  pctPool            Percentage of pooled environmental component\n"+"  msCoeffInHost            Parameter related to microbe selection strength\n"+"  msCoeffInEnv            Parameter related to microbe selection strength\n"+
                "  hsCoeff            Parameter related to host selection strength\n"+"  HMS_or_TMS         String HMS or TMS to specify host-mediated or trait-mediated microbe selection\n"
                + "\noptional arguments:\n";
        String footer = "\n";

        formatter.setWidth(80);


        CommandLineParser parser = new DefaultParser();
        CommandLine cmd = null;

        try {
            cmd = parser.parse(options, args);
            String[] pct_config = cmd.getArgs();

            if (cmd.hasOption("h") || args.length == 0) {
                formatter.printHelp(syntax, header, options, footer, true);
                System.exit(0);
            }
            if(cmd.hasOption("v")){
                System.out.println("Microbiosima "+VERSION);
                System.exit(0);
            }
            if (pct_config.length != 6){
                System.out.println("ERROR! Required exactly five argumennts for pct_env, pct_pool, msCoeff, hsCoeff and HMS_or_TMS. It got "
                        + pct_config.length + ": " + Arrays.toString(pct_config));
                formatter.printHelp(syntax, header, options, footer, true);
                System.exit(3);
            }
            else{
                pctEnv = Double.parseDouble(pct_config[0]);
                pctPool = Double.parseDouble(pct_config[1]);
                msCoeffInHost = Double.parseDouble(pct_config[2]);
                msCoeffInEnv = Double.parseDouble(pct_config[3]);
                hsCoeff = Double.parseDouble(pct_config[4]);
                if (pct_config[5].equals("HMS")) HMS_or_TMS=true;
                if (pct_config[5].equals("TMS")) HMS_or_TMS=false;
                if(pctEnv<0 || pctEnv >1){
                    System.out.println(
                            "ERROR: pctEnv (Percentage of environmental acquisition) must be between 0 and 1 (pctEnv="
                                    + pctEnv + ")! EXIT");
                    System.exit(3);
                }
                if(pctPool<0 || pctPool >1){
                    System.out.println(
                            "ERROR: pctPool (Percentage of pooled environmental component must) must be between 0 and 1 (pctPool="
                                    + pctPool + ")! EXIT");
                    System.exit(3);
                }
                if(msCoeffInHost<1){
                    System.out.println(
                            "ERROR: msCoeffInHost (parameter related to microbe selection strength) must be not less than 1 (msCoeff="
                                    + msCoeffInHost + ")! EXIT");
                    System.exit(3);
                }
                if(msCoeffInEnv<1){
                    System.out.println(
                            "ERROR: msCoeffInEnv (parameter related to microbe selection strength) must be not less than 1 (msCoeff="
                                    + msCoeffInEnv + ")! EXIT");
                    System.exit(3);
                }
                if(hsCoeff<1){
                    System.out.println(
                            "ERROR: hsCoeff (parameter related to host selection strength) must be not less than 1 (hsCoeff="
                                    + hsCoeff + ")! EXIT");
                    System.exit(3);
                }
                if (!(pct_config[5].equals("HMS")||pct_config[5].equals("TMS"))){
                    System.out.println(
                            "ERROR: HMS_or_TMS (parameter specifying host-mediated or trait-mediated selection) must be either 'HMS' or 'TMS' (HMS_or_TMS="
                                    + pct_config[5] + ")! EXIT");
                    System.exit(3);
                }

            }
            if (cmd.hasOption("config")){
                String[] configs = cmd.getOptionValues("config");
                populationSize = Integer.parseInt(configs[0]);
                microSize = Integer.parseInt(configs[1]);
                numberOfSpecies = Integer.parseInt(configs[2]);
                numberOfGeneration = Integer.parseInt(configs[3]);
                Ngene=Integer.parseInt(configs[4]);
                Ngenepm=Double.parseDouble(configs[5]);
                if (Ngenepm>Ngene){
                    System.out.println(
                            "ERROR: number of traits per microbe must not be greater than number of total traits! EXIT");
                    System.exit(3);
                }
            }
            if (cmd.hasOption("obs")){
                numberOfObservation= Integer.parseInt(cmd.getOptionValue("obs"));
            }
            if (cmd.hasOption("rep")){
                numberOfReplication= Integer.parseInt(cmd.getOptionValue("rep"));
            }



        } catch (ParseException e) {
            e.printStackTrace();
            System.exit(3);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("Configuration Summary:")
                .append("\n\tPopulation size: ").append(populationSize)
                .append("\n\tMicrobe size: ").append(microSize)
                .append("\n\tNumber of species: ").append(numberOfSpecies)
                .append("\n\tNumber of generation: ").append(numberOfGeneration)
                .append("\n\tNumber generation for observation: ").append(numberOfObservation)
                .append("\n\tNumber of replication: ").append(numberOfReplication)
                .append("\n\tNumber of total traits: ").append(Ngene)
                .append("\n\tNumber of traits per microbe: ").append(Ngenepm)
                .append("\n");
        System.out.println(sb.toString());



        double[] environment=new double[numberOfSpecies];//将环境中的微生物OTU占比设为等概率
        for (int i=0;i<numberOfSpecies;i++){
            environment[i]=1/(double)numberOfSpecies;
        }
        int[] fitnessToHost=new int[Ngene];//初始化25种性状
        int[] fitnessToMicrobeInHost=new int[Ngene];//初始化25种性状
        int[] fitnessToMicrobeInEnv=new int[Ngene];//初始化25种性状

        for (int rep=0;rep<numberOfReplication;rep++){
            String prefix = ""+(rep+1)+"_";
            String sufix;
            if (HMS_or_TMS)
                sufix = "_E" + pctEnv + "_P" + pctPool +"_HS"+hsCoeff+"_HMS"+ msCoeffInHost+"_EMS"+ msCoeffInEnv +".txt";
            else
                sufix = "_E" + pctEnv + "_P" + pctPool +"_HS"+hsCoeff+"_TMS"+ msCoeffInHost+"_EMS"+ msCoeffInEnv +".txt";
            System.out.println("Output 5 result files in the format of: "+prefix+"[****]" +sufix);
            try{
                PrintWriter file1= new PrintWriter(new BufferedWriter(new FileWriter(prefix+"gamma_diversity"+sufix)));
                PrintWriter file2= new PrintWriter(new BufferedWriter(new FileWriter(prefix+"alpha_diversity"+sufix)));
                PrintWriter file3= new PrintWriter(new BufferedWriter(new FileWriter(prefix+"beta_diversity"+sufix)));
                PrintWriter file4= new PrintWriter(new BufferedWriter(new FileWriter(prefix+"sum"+sufix)));
                PrintWriter file5= new PrintWriter(new BufferedWriter(new FileWriter(prefix+"inter_generation_distance"+sufix)));
                PrintWriter file6= new PrintWriter(new BufferedWriter(new FileWriter(prefix+"environment_population_distance"+sufix)));
                PrintWriter file7= new PrintWriter(new BufferedWriter(new FileWriter(prefix+"host_fitness"+sufix)));
                PrintWriter file8= new PrintWriter(new BufferedWriter(new FileWriter(prefix+"cos_theta"+sufix)));
                PrintWriter file9= new PrintWriter(new BufferedWriter(new FileWriter(prefix+"host_fitness_distribution"+sufix)));
                PrintWriter file10= new PrintWriter(new BufferedWriter(new FileWriter(prefix+"microbiome_fitness_in_host_distribution"+sufix)));
                PrintWriter file11= new PrintWriter(new BufferedWriter(new FileWriter(prefix+"bacteria_contents"+sufix)));
                PrintWriter file12= new PrintWriter(new BufferedWriter(new FileWriter(prefix+"individual_bacteria_contents"+sufix)));
                PrintWriter file13= new PrintWriter(new BufferedWriter(new FileWriter(prefix+"k"+sufix)));
                PrintWriter file14= new PrintWriter(new BufferedWriter(new FileWriter(prefix+"environmental_alpha_diversity"+sufix)));
                PrintWriter file15= new PrintWriter(new BufferedWriter(new FileWriter(prefix+"environmental_bacteria_contents_for_host"+sufix)));
                PrintWriter file16= new PrintWriter(new BufferedWriter(new FileWriter(prefix+"environmental_bacteria_contents_for_env"+sufix)));
                PrintWriter file17= new PrintWriter(new BufferedWriter(new FileWriter(prefix+"microbiome_fitness_in_environment_distribution"+sufix)));
                for (int i=0;i<Ngene;i++){
                    fitnessToMicrobeInHost[i]=MathUtil.getNextInt(2)-1;//随机分配微生物的25种性状数值-1、0、1.
                    fitnessToMicrobeInEnv[i]=MathUtil.getNextInt(2)-1;//随机分配微生物的25种性状数值-1、0、1.
                    fitnessToHost[i]=MathUtil.getNextInt(2)-1;//随机分配宿主的25种性状数值-1、0、1.
                }
                MathUtil.setSeed(rep%numberOfReplication);//生成随机数种子，防止每次的随机数都一样
                //SelectiveSpeciesRegistry根据微生物种类，性状总数，每个微生物性状总数，sm参数等赋予了每种OTU性状，并记录到geneCompositionRecords中
                SelectiveSpeciesRegistry ssr=new SelectiveSpeciesRegistry(numberOfSpecies,Ngene, Ngenepm, msCoeffInHost,msCoeffInEnv,fitnessToHost,fitnessToMicrobeInHost,fitnessToMicrobeInEnv);
                MathUtil.setSeed();//?
                SelectivePopulationE population=new SelectivePopulationE(microSize, environment, populationSize, pctEnv , pctPool,0,0,ssr,hsCoeff,HMS_or_TMS);

                while (population.getNumberOfGeneration()<numberOfGeneration){
                    population.sumSpecies();
                    if(population.getNumberOfGeneration()%numberOfObservation==0||population.getNumberOfGeneration()==1){
                        //file1.print(population.gammaDiversity(false));
                        //file2.print(population.alphaDiversity(false));
                        //file1.print("\t");
                        //file2.print("\t");
                        file1.println(population.gammaDiversity(true));
                        file2.println(population.alphaDiversity(true));
                        //file3.print(population.betaDiversity(true));
                        //file3.print("\t");
                        file3.println(population.BrayCurtis(true));
                        file4.println(population.printOut());//微生物distribution，每种OTU的占比
                        file5.println(population.interGenerationDistance());
                        file6.println(population.environmentPopulationDistance());
                        file7.print(population.averageHostFitness());//宿主平均适应度
                        file7.print("\t");
                        file7.println(population.varianceHostFitness());//宿主适应度的方差
                        file8.println(population.cosOfMH());
                        file9.println(population.printOutHFitness());//每种otu对宿主的fitness
                        file10.println(population.printOutMFitnessInHost());//每种otu自身的fitness
                        file17.println(population.printOutMFitnessInEnvironment());//每种otu自身的fitness
                        file11.println(population.printBacteriaContents());
                        file13.println(population.printk());
                        file14.println(population.alphaDiversityInEnvironment());
                        file15.println(population.printBacteriaContentsInEnvForHost());
                        file16.println(population.printBacteriaContentsInEnv());
                    }
                    population.getNextGen();
                }
                for (SelectiveIndividualE host:population.getIndividuals()){
                    file12.println(host.printBacteriaContents());
                }
                file1.close();
                file2.close();
                file3.close();
                file4.close();
                file5.close();
                file6.close();
                file7.close();
                file8.close();
                file9.close();
                file10.close();
                file11.close();
                file12.close();
                file13.close();
                file14.close();
                file15.close();
                file16.close();
                file17.close();

            }catch (IOException e) {
                e.printStackTrace();
            }
        }



    }

}


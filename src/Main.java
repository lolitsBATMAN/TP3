/**
 * Chengzong Jiang (20122046)
 * Ahmad El-Rahim (20157008)
 */

import java.io.*;
import java.util.Scanner;
import java.util.ArrayList;

public class Main {
    // 1.1
    public static ArrayList<DBentry> readFasta(String file) throws FileNotFoundException {
        ArrayList<DBentry> database =  new ArrayList<>();
        try (Scanner sc = new Scanner(new File(file))) {

            while (sc.hasNextLine()) {
                String infos = sc.nextLine();
                String sequence = sc.nextLine();

                //create entry that contains these 2 strings
                DBentry currEntry = new DBentry(infos,sequence);

                //add entry to database
                database.add(currEntry);
            }
        }
        return database;
    }

    //1.2 and 1.3
    //la position de chaque KMER est l'index dans le tableau
    public static ArrayList<String> kmer (String word, int k){
        ArrayList<String> mer = new ArrayList<String>();

        for (int i  = 0; i <= word.length() - k ; i++){
            mer.add(word.substring(i, i+k));
        }

        return mer;
    }

    public static ArrayList<ArrayList<String>> hsp(String input, ArrayList<DBentry> data, String seed, ArrayList<ArrayList<Integer>> position, ArrayList<ArrayList<Integer>> positioninput){
        int k = seed.length();
        ArrayList<ArrayList<String>> hsps = new ArrayList<>();
        ArrayList<String> kinput = kmer(input, k);

        for (int y = 0; y < data.size(); y++) {
            ArrayList<String> kdata = kmer(data.get(y).getSequence(), k);
            hsps.add(new ArrayList<>());
            position.add(new ArrayList<>());
            positioninput.add(new ArrayList<>());
            for (int i = 0; i < kinput.size(); i++) {
                for (int j = 0; j < kdata.size(); j++) {
                    boolean bruh = true;
                    for (int z = 0; z < k; z++) {
                        if (kinput.get(i).charAt(z) != kdata.get(j).charAt(z) && seed.charAt(z) == '1') {
                            bruh = false;
                            break;
                        }
                    }
                    if (bruh) {
                        hsps.get(y).add(kinput.get(i));
                        position.get(y).add(j); // y ieme sequence dans la database , j ieme kmer = index de cette sequence
                        positioninput.get(y).add(i);// y ieme sequence dans la database, i eme kmer = index de notre sequence inconnue
                    }
                }
            }
        }
        return hsps;
    }

    public static void removeDuplicates(ArrayList<ArrayList<String>> list, ArrayList<ArrayList<Integer>> positionInput, ArrayList<ArrayList<Integer>> positionDatabase) {

        // Create a new ArrayList
        ArrayList<ArrayList<String>> newList = new ArrayList<>();
        ArrayList<ArrayList<Integer>> positionIn = new ArrayList<>();
        ArrayList<ArrayList<Integer>> positionDb = new ArrayList<>();

        // Traverse through the first list
        for (int i = 0;  i < list.size(); i++) {
            ArrayList<String> newist = new ArrayList<>();
            ArrayList<Integer> posIn = new ArrayList<>();
            ArrayList<Integer> posDb = new ArrayList<>();
            if (!list.get(i).isEmpty()) {
                for (int j = 0; j < list.get(i).size(); j++) {

                    if (!newist.contains(list.get(i).get(j)) || !posIn.contains(positionInput.get(i).get(j)) || !posDb.contains(positionDatabase.get(i).get(j))) {
                        newist.add(list.get(i).get(j));
                        posIn.add(positionInput.get(i).get(j));
                        posDb.add(positionDatabase.get(i).get(j));
                    }
                }
            }
            newList.add(newist);
            positionIn.add(posIn);
            positionDb.add(posDb);
        }

        list.clear();
        positionInput.clear();
        positionDatabase.clear();

        list.addAll(newList);
        positionInput.addAll(positionIn);
        positionDatabase.addAll(positionDb);
    }

    //1.4 extend des hsps
    public static ArrayList<ArrayList<String>> glouton (ArrayList<ArrayList<String>> hsp, ArrayList<DBentry> database, String input, ArrayList<ArrayList<Integer>> positionInput, ArrayList<ArrayList<Integer>> positionDatabase, Integer seuil, String seed){

        for (int i = 0; i < hsp.size(); i++){
            if (!hsp.get(i).isEmpty()) {
                for (int j = 0; j < hsp.get(i).size(); j++) {

                    int seqInputFirst = positionInput.get(i).get(j)-1;
                    int seqInputLast = positionInput.get(i).get(j) + seed.length();
                    int seqDataFirst = positionDatabase.get(i).get(j)-1;
                    int seqDataLast = positionDatabase.get(i).get(j) + seed.length();

                    int maxScore = 0;
                    int score = 0;

                    boolean extendLeft = true;
                    boolean extendRight = true;

                    while(extendLeft || extendRight){
                        int scoreLeft = 0;
                        int scoreRight = 0;

                        if (seqInputFirst == -1 || seqDataFirst == -1){
                            extendLeft = false;
                        } else {
                            if (input.charAt(seqInputFirst) == database.get(i).getSequence().charAt(seqDataFirst)){
                                scoreLeft += 5;
                            } else {
                                scoreLeft -= 4;
                            }
                        }

                        if (seqDataLast == database.get(i).getSequence().length()   || seqInputLast == input.length()) {
                            extendRight = false;
                        } else {
                            if (input.charAt(seqInputLast) == database.get(i).getSequence().charAt(seqDataLast)){
                                scoreRight += 5;
                            } else {
                                scoreRight -= 4;
                            }
                        }

                        if (!extendLeft && extendRight) {
                            score += scoreRight;
                            if (maxScore - score < seuil){
                                hsp.get(i).set(j, hsp.get(i).get(j)+input.charAt(seqInputLast));
                                seqInputLast++;
                                seqDataLast++;
                            }
                        }

                        if (extendLeft && !extendRight) {
                            score += scoreLeft;
                            if (maxScore - score < seuil) {
                                hsp.get(i).set(j, input.charAt(seqInputFirst) + hsp.get(i).get(j));
                                positionDatabase.get(i).set(j, positionDatabase.get(i).get(j)-1);
                                positionInput.get(i).set(j, positionInput.get(i).get(j)-1);
                                seqInputFirst--;
                                seqDataFirst--;
                            }
                        }

                        if (extendLeft && extendRight) {
                            score += scoreLeft;
                            if (scoreLeft >= scoreRight){
                                if (maxScore - score < seuil) {
                                    hsp.get(i).set(j, input.charAt(seqInputFirst) + hsp.get(i).get(j));
                                    positionDatabase.get(i).set(j, positionDatabase.get(i).get(j)-1);
                                    positionInput.get(i).set(j, positionInput.get(i).get(j)-1);
                                    seqInputFirst--;
                                    seqDataFirst--;
                                }
                            } else {
                                score += scoreRight;
                                if (maxScore - score < seuil) {
                                    hsp.get(i).set(j, hsp.get(i).get(j) + input.charAt(seqInputLast));
                                    seqInputLast++;
                                    seqDataLast++;
                                }
                            }
                        }

                        if (score > maxScore){
                            maxScore = score;
                        }

                        if (maxScore - score >= seuil){
                            extendLeft = false;
                            extendRight = false;
                        }
                    }
                }
            }
        }

        removeDuplicates(hsp, positionInput, positionDatabase);

        return hsp;
    }
    //1.5
    public static ArrayList<ArrayList<String>> fusionHsp (String input, ArrayList<ArrayList<String>> extendedHsps, ArrayList<ArrayList<Integer>> positionInput, ArrayList<ArrayList<Integer>> positionDatabase) {

        ArrayList<ArrayList<String>> result = new ArrayList<>();
        ArrayList<ArrayList<Integer>> in = new ArrayList<>();
        ArrayList<ArrayList<Integer>> db = new ArrayList<>();

        for (int i = 0; i < positionDatabase.size(); i++) {

            ArrayList<Integer> offset = new ArrayList<>();
            ArrayList<ArrayList<String>> hsps = new ArrayList<>();
            ArrayList<ArrayList<Integer>> positionIn = new ArrayList<>();
            ArrayList<ArrayList<Integer>> positionDb = new ArrayList<>();

            if (!extendedHsps.get(i).isEmpty()){
                for (int j = 0; j < positionDatabase.get(i).size(); j++) {

                    int nb = positionInput.get(i).get(j) - positionDatabase.get(i).get(j);

                    if (!offset.contains(nb)) {
                        offset.add(nb);
                        hsps.add(new ArrayList<>());
                        hsps.get(hsps.size()-1).add(extendedHsps.get(i).get(j));
                        positionIn.add(new ArrayList<>());
                        positionIn.get(hsps.size()-1).add(positionInput.get(i).get(j));
                        positionDb.add(new ArrayList<>());
                        positionDb.get(hsps.size()-1).add(positionDatabase.get(i).get(j));
                    } else {
                        hsps.get(offset.indexOf(nb)).add(extendedHsps.get(i).get(j));
                        positionIn.get(offset.indexOf(nb)).add(positionInput.get(i).get(j));
                        positionDb.get(offset.indexOf(nb)).add(positionDatabase.get(i).get(j));
                    }
                }

                for (int k = 0; k < hsps.size(); k++){
                    for (int l = 0; l < hsps.get(k).size(); l++){
                        for (int m = 0; m < hsps.get(k).size(); m++){
                            if (m != l){
                                int positionFirstInit = positionIn.get(k).get(l);
                                int positionFirstFinal = positionFirstInit + hsps.get(k).get(l).length() - 1;
                                int positionNextInit = positionIn.get(k).get(m);
                                int positionNextFinal = positionNextInit + hsps.get(k).get(m).length() - 1;

                                if (positionNextInit  <= positionFirstFinal && positionFirstFinal <= positionNextFinal){
                                    hsps.get(k).set(l, input.substring(positionFirstInit, positionNextFinal));
                                    hsps.get(k).remove(m);
                                    positionIn.get(k).remove(m);
                                    positionDb.get(k).remove(m);
                                }
                            }
                        }
                    }
                }
            }
            result.add(new ArrayList<>());
            in.add(new ArrayList<>());
            db.add(new ArrayList<>());
            for (int k = 0; k < hsps.size(); k++){
                for (int l = 0; l < hsps.get(k).size(); l++){
                    result.get(i).add(hsps.get(k).get(l));
                    in.get(i).add(positionIn.get(k).get(l));
                    db.get(i).add(positionDb.get(k).get(l));
                }
            }
        }
        positionDatabase.clear();
        positionDatabase.addAll(db);
        positionInput.clear();
        positionInput.addAll(in);

        return result;
    }
    //1.6
    public static int databaseTotalSize(ArrayList<DBentry> database){
        int len = 0;
        for (DBentry dBentry : database) {
            len += dBentry.sequence.length();
        }
        return len;
    }

    public static void calcul(ArrayList<ArrayList<HspStat>> resultat, ArrayList<HspStat> maximum,String input, ArrayList<DBentry> database,ArrayList<ArrayList<String>> hsps, ArrayList<ArrayList<Integer>> positionInput, ArrayList<ArrayList<Integer>> positionDatabase, double seuil2){
        ArrayList<ArrayList<HspStat>> result = new ArrayList<>();

        int m = databaseTotalSize(database);
        for(int i = 0; i < positionDatabase.size(); i++){
            result.add(new ArrayList<>());
            if (!hsps.get(i).isEmpty()){
                for (int j = 0 ;j < positionDatabase.get(i).size(); j++){
                    result.get(i).add(new HspStat(input,hsps.get(i).get(j), positionInput.get(i).get(j),positionDatabase.get(i).get(j), database.get(i), m));
                }
            }
        }

        ArrayList<HspStat> maxi = new ArrayList<>();

        for (ArrayList<HspStat> hspStats : result) {
            double max = -1;
            HspStat maxim = null;
            for (int j = 0; j < hspStats.size(); j++) {
                if (hspStats.get(j).geteValue() >= seuil2) {
                    hspStats.remove(j);
                } else {
                    if (hspStats.get(j).getBitScore() >= max) {
                        max = hspStats.get(j).getBitScore();
                        maxim = hspStats.get(j);
                    }
                }
            }
            maxi.add(maxim);
        }
        resultat.addAll(result);
        maximum.addAll(maxi);
    }


    //1.7
    public static void plast(String input, String output, int e, double ss, String seed) throws FileNotFoundException {
        ArrayList<DBentry> database = readFasta(output);
        ArrayList<ArrayList<Integer>> position = new ArrayList<>();
        ArrayList<ArrayList<Integer>> positionInput = new ArrayList<>();
        ArrayList<ArrayList<String>> hsps = hsp(input, database,seed,position,positionInput);
        ArrayList<ArrayList<String>> hspExtend = glouton(hsps, database, input, positionInput, position, e, seed);
        ArrayList<ArrayList<String>> fuse = fusionHsp (input, hspExtend, positionInput, position);
        ArrayList<HspStat> max = new ArrayList<>();
        ArrayList<ArrayList<HspStat>> result = new ArrayList<>();
        calcul(result, max, input, database, fuse, positionInput,position, ss);
        ArrayList<HspStat> sort = new ArrayList<>();

        for (HspStat hspStat : max) {
            if (hspStat != null) {
                sort.add(hspStat);
            }
        }

        boolean sorted = false;
        HspStat temp;
        while(!sorted){
            sorted = true;
            for (int i = 0; i < sort.size() - 1; i++){
                if (sort.get(i).getBitScore() < sort.get(i+1).getBitScore()){
                    temp = sort.get(i);
                    sort.set(i, sort.get(i+1));
                    sort.set(i+1, temp);
                    sorted = false;
                }
            }
        }

        if (!sort.isEmpty()){
            System.out.println();
            System.out.println("---------------------------------------------");
            System.out.println("Input: "+input);
            System.out.println();
            for (HspStat i : sort){
                System.out.println(i.getDatabase().getSequenceInfo());
                System.out.println("# Best HSP score: " + i.getBruteScore() + ", bitscore: "+i.getBitScore()+ ", evalue: "+ i.geteValue());
                System.out.println(i.getPositionInput()+ " "+i.getCurrentHsp() + " "+ (i.getPositionInput()+i.getCurrentHsp().length()-1));
                System.out.println(i.getPositionDatabase() + " "+i.getDatabase().getSequence().substring(i.getPositionDatabase(), i.getPositionDatabase()+i.getCurrentHsp().length())+" "+ (i.getPositionDatabase()+i.getCurrentHsp().length()-1));
                System.out.println();
            }
            System.out.println("---------------------------------------------");
            System.out.println("total: " + sort.size());
        } else {
            System.out.println("No results");
        }

    }

    public static void main(String[] args) throws IOException {

        BufferedReader in = new BufferedReader(new InputStreamReader(System.in));

        System.out.println("Entrer le input: ");

        String input = in.readLine();

        System.out.println("Entrer le path du database: ");

        String database = in.readLine();

        System.out.println("Entrer la seuil E: ");

        int e = Integer.parseInt(in.readLine());

        System.out.println("Entrer la seuil SS: ");

        double ss = Double.parseDouble(in.readLine());

        System.out.println("Entrer la graine: ");

        String seed = in.readLine();

        plast(input, database, e, ss, seed);
        //plast("src/unknown.fasta", "src/tRNAs.fasta", 4, 0.001,"111010010100110111");
    }
}

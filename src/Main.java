import javax.lang.model.type.NullType;
import java.io.*;
import java.util.Scanner;
import java.util.ArrayList;

public class Main {
    // 1.1
    //je m'en criss un peu ici so j'ai tout crisser dans un tableau, dans le fond
    //les index 0,2,4,6... sont les >dsad|dasd|sdads| et les index 1,3,5,... sont
    //les nucleotides, good luck
    //trop de trouble de les filtrer

    //2e argument la dabatase quon creer dans la main
    public static ArrayList<DBentry> readFasta(String file) throws FileNotFoundException {
        ArrayList<DBentry> database =  new ArrayList<DBentry>();
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

    //retourne une list avec les hsps
    public static ArrayList<ArrayList<String>> hsp(String input, ArrayList<DBentry> data, String seed, ArrayList<ArrayList<Integer>> position, ArrayList<ArrayList<Integer>> positioninput){
        int k = seed.length();
        ArrayList<ArrayList<String>> hsps = new ArrayList<ArrayList<String>>();
        ArrayList<String> kinput = kmer(input, k);

        for (int y = 0; y < data.size(); y++) {
            ArrayList<String> kdata = kmer(data.get(y).getSequence(), k);
            hsps.add(new ArrayList<String>());
            position.add(new ArrayList<Integer>());
            positioninput.add(new ArrayList<Integer>());
            for (int i = 0; i < kinput.size(); i++) {
                for (int j = 0; j < kdata.size(); j++) {
                    boolean bruh = true;
                    for (int z = 0; z < k; z++) {
                        if (kinput.get(i).charAt(z) != kdata.get(j).charAt(z) && seed.charAt(z) == '1') {
                            bruh = false;
                            break;
                        } else if (kinput.get(i).charAt(z) == kdata.get(j).charAt(z) && seed.charAt(z) == '0') {
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
        ArrayList<ArrayList<String>> newList = new ArrayList<ArrayList<String>>();
        ArrayList<ArrayList<Integer>> positionIn = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<Integer>> positionDb = new ArrayList<ArrayList<Integer>>();

        // Traverse through the first list
        for (int i = 0;  i < list.size(); i++) {
            ArrayList<String> newist = new ArrayList<String>();
            ArrayList<Integer> posIn = new ArrayList<Integer>();
            ArrayList<Integer> posDb = new ArrayList<Integer>();
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

    //1.4
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


//        System.out.println("INITAL SIZE");
//        System.out.println(positionInput.get(0).size());

        ArrayList<ArrayList<String>> result = new ArrayList<>();
        ArrayList<ArrayList<Integer>> in = new ArrayList<>();
        ArrayList<ArrayList<Integer>> db = new ArrayList<>();

        for (int i = 0; i < positionDatabase.size(); i++) {
            ArrayList<Integer> offset = new ArrayList<>();
            ArrayList<ArrayList<String>> hsps = new ArrayList<>();
            ArrayList<ArrayList<Integer>> positionIn = new ArrayList<>();
            ArrayList<ArrayList<Integer>> positionDb = new ArrayList<>();

//            System.out.println(positionInput);
//            System.out.println(positionDatabase);

            //System.out.println(!extendedHsps.get(i).isEmpty());
            if (!extendedHsps.get(i).isEmpty()){
//                System.out.println(positionInput);
//                System.out.println(positionDatabase);
                for (int j = 0; j < positionDatabase.get(i).size(); j++) {
//                    System.out.println(positionInput);
//                    System.out.println(positionDatabase);
//                    System.out.println();
//                    System.out.println(extendedHsps.get(i));
//                    System.out.println(positionDatabase.get(i));
//                    System.out.println(positionInput.get(i));


                    int nb = positionInput.get(i).get(j) - positionDatabase.get(i).get(j);
                    //System.out.println(nb);
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

//                if (i == 0){
//                    System.out.println("OFFSET");
//                    System.out.println(offset);
//                    System.out.println(offset.size());
//                }

                for (int k = 0; k < hsps.size(); k++){
                    for (int l = 0; l < hsps.get(k).size(); l++){
                        for (int m = 0; m < hsps.get(k).size(); m++){
                            if (m != l){
                                int positionFirstInit = positionIn.get(k).get(l);
                                int positionFirstFinal = positionFirstInit + hsps.get(k).get(l).length() - 1;
                                int positionNextInit = positionIn.get(k).get(m);
                                int positionNextFinal = positionNextInit + hsps.get(k).get(m).length() - 1;

                                if (positionNextInit  <= positionFirstFinal && positionFirstFinal <= positionNextFinal){
                                    //System.out.println(input.length());
                                    //System.out.println(input);
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

//        System.out.println("FINAL SIZE");
//        System.out.println(positionInput.get(0).size());

        return result;
    }
    //1.6
    public static int databaseTotalSize(ArrayList<DBentry> database){
        int len = 0;
        for (int i= 0; i<database.size(); i++){
            len += database.get(i).sequence.length();
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
                    //System.out.println(positionDatabase.get(i).size());
                    result.get(i).add(new HspStat(input,hsps.get(i).get(j), positionInput.get(i).get(j),positionDatabase.get(i).get(j), database.get(i).getSequence(), m));
                }
            }
        }

        ArrayList<HspStat> maxi = new ArrayList<>();

        for (int i = 0; i < result.size(); i++){
            double max = -1;
            HspStat maxim = null;
            for (int j = 0; j < result.get(i).size(); j++){
                System.out.println(result.get(i).get(j).geteValue());
                System.out.println(result.get(i).get(j).getBruteScore());
                //System.out.println(result.get(i).get(j).getBitScore());
                //System.out.println(result.get(i).get(j));
                if (result.get(i).get(j).geteValue() >= seuil2){
                    result.get(i).remove(j);
                } else {
                    if (result.get(i).get(j).getBitScore() >= max){
                        max = result.get(i).get(j).getBitScore();
                        maxim = result.get(i).get(j);
                    }
                }
            }
            maxi.add(maxim);
        }
        resultat.addAll(result);
        maximum.addAll(maxi);
    }


    //1.7
    public static void plast(String input, String output, String seed, int seuil, double seuil2) throws FileNotFoundException {
        ArrayList<DBentry> unknown = readFasta(input);
        ArrayList<DBentry> database = readFasta(output);
        ArrayList<ArrayList<Integer>> position = new ArrayList<>();
        ArrayList<ArrayList<Integer>> positionInput = new ArrayList<>();
        ArrayList<ArrayList<String>> hsps = hsp(unknown.get(0).getSequence(), database,seed,position,positionInput);
        ArrayList<ArrayList<String>> hspExtend = glouton(hsps, database, unknown.get(0).getSequence(), positionInput, position, seuil, seed);
        ArrayList<ArrayList<String>> fuse = fusionHsp (unknown.get(0).getSequence(), hspExtend, positionInput, position);
        System.out.println(position);
        ArrayList<HspStat> max = new ArrayList<>();
        ArrayList<ArrayList<HspStat>> result = new ArrayList<>();
        calcul(result, max, unknown.get(0).getSequence(), database, fuse, positionInput,position, seuil2);

        System.out.println(result);
        System.out.println(max);
    }

    public static void main(String[] args) throws FileNotFoundException {
//        ArrayList<DBentry> db = new ArrayList<DBentry>();
//        DBentry x = new DBentry("BRUH","ABCXXXXABC");
//        db.add(x);
//
//        String input = "PEFABCEFP";
//        ArrayList<ArrayList<Integer>> position = new ArrayList<ArrayList<Integer>>();
//        ArrayList<ArrayList<Integer>> positionInput = new ArrayList<ArrayList<Integer>>();
//        ArrayList<ArrayList<String>> hsps = hsp(input, db,"111",position,positionInput);
//        System.out.println(glouton(hsps, db,input,positionInput, position, 4,"111"));
//        System.out.println(hsps);
//        System.out.println(positionInput);
//        System.out.println(position);

        plast("src/unknown.fasta","src/tRNAs.fasta", "11111111111", 5, 10);
    }
}

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace OCCFRecSys
{
    public class BalanceRWR_Recommender : Recommender
    {
        public int fold;
        public double crit;
        string candidate_item;
        // RWR variables
        int num_neighbor;                                               // graph 생성 시 파라미터 (이웃 수)
        double  e, damping_factor;
        public int user_size;                                           // graph의 user 수
        List<WeightedLink> links_like, links_dislike;                   // graph의 links
        public Dictionary<int, List<Score>> recommend_user_item;        // 타겟 유저들에게 추천할 아이템들과 점수 저장     

        public BalanceRWR_Recommender(string candidate_item, int num_recommend, int fold, double crit, StreamReader training, StreamReader test)
        {
            this.candidate_item = candidate_item;
            this.fold = fold;
            this.crit = crit;

            Init(num_recommend, training, test);
        }

        public override void Init(int num_recommend, StreamReader training, StreamReader test)
        {
            num_neighbor = 50;  // graph 생성 시 파라미터 (이웃 수)
            // RWR variables
            e = 0.0d;
            damping_factor = 0.85d;
            this.num_recommend = num_recommend;

            // test set 읽기            
            ReadTestSet(test, candidate_item);
            test.Close();
            Console.WriteLine("Read test set: Complete");
            
            // training set 읽기
            ReadTrainingSet(training);
            training.Close();
            Console.WriteLine("Read training set: Complete");
            
        }

        public override void ReadTrainingSet(StreamReader sr)
        {
            links_like = new List<WeightedLink>();
            links_dislike = new List<WeightedLink>();

            while (!sr.EndOfStream)
            {
                line = sr.ReadLine().Split(new char[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                int user_id = int.Parse(line[0].ToString());
                int item_id = int.Parse(line[1].ToString()) + Program.num_user;
                int rating = int.Parse(line[2].ToString());

                if (rating == 1)
                    links_like.Add(new WeightedLink(user_id, item_id, 1));
                else
                    links_dislike.Add(new WeightedLink(user_id, item_id, 1));
            }
        }

        public override void Recommend()
        {
            // RWR 점수 저장
            recommend_user_item = new Dictionary<int, List<Score>>();
            
            WeightedUndirectedGraph positive_graph = new WeightedUndirectedGraph(links_like, num_neighbor, e);
            WeightedUndirectedGraph negative_graph = new WeightedUndirectedGraph(links_dislike, num_neighbor, e);
            user_size = positive_graph.user_size;

            for (int i = 1; i < user_size; i++)
            {
                // 해당 사용자가 test user가 아니면 스킵
                if (!correct_user_item.ContainsKey(i)) continue;
                else if (i % 50 == 0)
                    Console.WriteLine("{0} Recommendation: user {1} / {2}", fold, i, (user_size - 1));
                
                Dictionary<int, double> like_items = PerformRWR(positive_graph, i);
                Dictionary<int, double> dislike_items = PerformRWR(negative_graph, i);

                double dislike_min = 1d;
                foreach (int item_id in dislike_items.Keys)
                {
                    if (dislike_items[item_id] == 0) continue;
                    if (dislike_min > dislike_items[item_id])
                        dislike_min = dislike_items[item_id];
                }

                List<Score> recommend_items = new List<Score>();
                foreach (int item_id in like_items.Keys)
                {
                    double final = 0d;
                    //if(Program.remove_items != null)
                    //    if(Program.remove_items.ContainsKey(i))
                    //        if (Program.remove_items[i].Contains(item_id)) continue;

                    if (!dislike_items.ContainsKey(item_id)) continue;
                    final = like_items[item_id] - dislike_min;
                    //else if (dislike_items[item_id] == 0)
                    //{
                    //    final = like_items[item_id] - dislike_items[item_id]; // dislike_min;
                    //    //final = like_items[item_id] * (1 - dislike_min);
                    //}
                    //else
                    //{
                    //    final = like_items[item_id] - dislike_items[item_id];
                    //    //final = like_items[item_id] * (1 - dislike_items[item_id]);
                    //}

                    Score user_score = new Score(item_id, final);
                    recommend_items.Add(user_score);
                }
                recommend_items.Sort(Score.CompareScore);
                recommend_user_item.Add(i, recommend_items);
                
            }
            Console.WriteLine("Compute RWR scores of {0}: Complete", fold);
        }

        /// <summary>
        /// user-item 형태로 모델링된 bipartite graph를 사용하여 모든 test user에 대해 RWR 수행 
        /// </summary>
        /// <param name="model_graph">user-item에 대한 bipartite graph</param>
        /// <param name="sw">추천 결과(사용자 id, 아이템 id, 점수)가 출력될 파일</param>
        public Dictionary<int, double> PerformRWR(WeightedUndirectedGraph model_graph, int target)
        {
            HashSet<int> i_items = new HashSet<int>();
            if (model_graph.graph.ContainsKey(target))
            {
                HashSet<Weight> i_weight = model_graph.graph[target];
                foreach (Weight items in i_weight)
                    i_items.Add(items.id);
            }

            // RWR 수행
            RWR ranker = new RWR(model_graph, damping_factor, target);
            List<Weight> score = ranker.Calculate();

            double sum = 0d;            
            // RWR 점수를 기준으로 타겟 사용자가 선호할 만한 TOP N개의 아이템 출력
            Dictionary<int, double> recommend_items = new Dictionary<int, double>();
            foreach (Weight s in score)
                if (s.id >= user_size)
                    sum += s.w;

            foreach (Weight s in score)
            {
                if (s.id < user_size ) continue; // || i_items.Contains(s.id)

                if (candidate_item == "Longtail_items")
                    if (Program.longtails.Contains(s.id)) continue;

                recommend_items.Add(s.id, s.w / sum);
            }
            return recommend_items;
        }        

        public override void PrintResults(StreamWriter sw)
        {
            ComputeMRRatN(recommend_user_item, sw);
            //ComputeAUC(recommend_user_item, sw);
            ComputeHLU(recommend_user_item, sw);
            ComputeAccuary(recommend_user_item, sw);
            ComputeNDCG(recommend_user_item, sw);
            Console.WriteLine("Compute accuracy of recommendation: Complete");
        }
    }
}

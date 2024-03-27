#include "copyright.h"
#include "generalized_born.h"

namespace stormm {
namespace generalized_born_defaults {

//-------------------------------------------------------------------------------------------------
NeckGeneralizedBornTable::NeckGeneralizedBornTable() :
    table_size{21}, neck_cut{default_gb_neck_cut}, kscale{default_gb_kscale},
    neck_max_separation{static_cast<size_t>(table_size * table_size), "gb_table_maxsep"},
    neck_max_value{static_cast<size_t>(table_size * table_size), "gb_table_maxval"},
    sp_neck_max_separation{static_cast<size_t>(table_size * table_size), "gb_table_maxsep"},
    sp_neck_max_value{static_cast<size_t>(table_size * table_size), "gb_table_maxval"}
{
  // Load the maximum separations table for neck GB
  std::vector<double> sep_buffer = {
    2.26684999, 2.32547998, 2.38397002, 2.44234991, 2.50057006, 2.55867004, 2.61663008, 2.67443991,
    2.73212004, 2.78964996, 2.84704995, 2.90429997, 2.96141005, 3.01839995, 3.07523990, 3.13195992,
    3.18853998, 3.24498010, 3.30132008, 3.35752010, 3.41359997,
    2.31190991, 2.37017012, 2.42829990, 2.48632002, 2.54419994, 2.60196996, 2.65961003, 2.71710992,
    2.77449012, 2.83174992, 2.88887000, 2.94585991, 3.00272989, 3.05947995, 3.11610007, 3.17260003,
    3.22897005, 3.28521991, 3.34136009, 3.39738011, 3.45072007,
    2.35758996, 2.41548991, 2.47328997, 2.53097010, 2.58854008, 2.64599991, 2.70333004, 2.76056004,
    2.81766009, 2.87465000, 2.93151999, 2.98827004, 3.04489994, 3.10141993, 3.15781999, 3.21410990,
    3.27027988, 3.32633996, 3.38229990, 3.43812990, 3.49387002,
    2.40380001, 2.46138000, 2.51885009, 2.57623005, 2.63351011, 2.69067001, 2.74773002, 2.80468988,
    2.86152005, 2.91826010, 2.97488999, 3.03139997, 3.08781004, 3.14409995, 3.20030999, 3.25638008,
    3.31237006, 3.36824989, 3.42402005, 3.47970009, 3.53526998,
    2.45044994, 2.50773001, 2.56491995, 2.62200999, 2.67899990, 2.73589993, 2.79270005, 2.84940004,
    2.90598989, 2.96250010, 3.01889992, 3.07518005, 3.13138008, 3.18747997, 3.24346995, 3.29937005,
    3.35514998, 3.41085005, 3.46645999, 3.52196002, 3.57737994,
    2.49749994, 2.55450010, 2.61142993, 2.66825008, 2.72498989, 2.78163004, 2.83818007, 2.89463997,
    2.95100999, 3.00728989, 3.06346011, 3.11953998, 3.17553997, 3.23143005, 3.28723001, 3.34294009,
    3.39856005, 3.45409012, 3.50952005, 3.56487989, 3.62014008,
    2.54488993, 2.60163999, 2.65829992, 2.71487999, 2.77133989, 2.82780004, 2.88411999, 2.94034004,
    2.99650002, 3.05256009, 3.10853004, 3.16441989, 3.22021008, 3.27591991, 3.33154011, 3.38706994,
    3.44252992, 3.49789000, 3.55315995, 3.60836005, 3.66348004,
    2.59259009, 2.64910007, 2.70552993, 2.76187992, 2.81815004, 2.87434006, 2.93043995, 2.98645997,
    3.04240990, 3.09826994, 3.15404010, 3.20973992, 3.26536012, 3.32088995, 3.37632990, 3.43169999,
    3.48698997, 3.54219007, 3.59731007, 3.65236998, 3.70734000,
    2.64053988, 2.69684005, 2.75305009, 2.80918002, 2.86523008, 2.92122006, 2.97711992, 3.03294992,
    3.08870006, 3.14437008, 3.19995999, 3.25548005, 3.31090999, 3.36627007, 3.42156005, 3.47676992,
    3.53189993, 3.58695006, 3.64193010, 3.69684005, 3.75166988,
    2.68873000, 2.74482012, 2.80082989, 2.85676003, 2.91262007, 2.96841002, 3.02412009, 3.07976007,
    3.13532996, 3.19081998, 3.24622989, 3.30156994, 3.35684991, 3.41205001, 3.46718001, 3.52222991,
    3.57720995, 3.63212991, 3.68695998, 3.74173999, 3.79643989,
    2.73712993, 2.79302001, 2.84884000, 2.90458989, 2.96026993, 3.01587009, 3.07139993, 3.12685990,
    3.18225002, 3.23757005, 3.29281998, 3.34801006, 3.40313005, 3.45814991, 3.51314998, 3.56804991,
    3.62290001, 3.67767000, 3.73236990, 3.78700995, 3.84158993,
    2.78572011, 2.84142995, 2.89706993, 2.95264006, 3.00813007, 3.06356001, 3.11892009, 3.17422009,
    3.22946000, 3.28462005, 3.33971000, 3.39474010, 3.44970989, 3.50460005, 3.55943990, 3.61420989,
    3.66891003, 3.72356009, 3.77814007, 3.83263993, 3.88709998,
    2.83446002, 2.89000010, 2.94547009, 3.00088000, 3.05621004, 3.11146998, 3.16669011, 3.22182989,
    3.27689004, 3.33190989, 3.38685012, 3.44174004, 3.49656010, 3.55132008, 3.60601997, 3.66066003,
    3.71522999, 3.76975012, 3.82420993, 3.87859988, 3.93292999,
    2.88334990, 2.93873000, 2.99404001, 3.04928994, 3.10447001, 3.15959001, 3.21463990, 3.26962996,
    3.32455993, 3.37943006, 3.43424010, 3.48898005, 3.54365993, 3.59829998, 3.65286994, 3.70737004,
    3.76183009, 3.81622005, 3.87055993, 3.92483997, 3.97904992,
    2.93233991, 2.98760009, 3.04276991, 3.09786010, 3.15290999, 3.20787001, 3.26277995, 3.31764007,
    3.37242007, 3.42716002, 3.48183990, 3.53661990, 3.59100008, 3.64550996, 3.69994998, 3.75434995,
    3.80867004, 3.86295009, 3.91718006, 3.97133994, 4.02545023,
    2.98150992, 3.03660011, 3.09162998, 3.14658999, 3.20148993, 3.25632000, 3.31110001, 3.36580992,
    3.42047000, 3.47507000, 3.52962995, 3.58411002, 3.63855004, 3.69292998, 3.74725008, 3.80152988,
    3.85575008, 3.90990996, 3.96403003, 4.01808977, 4.07211018,
    3.03074002, 3.08571005, 3.14060998, 3.19543004, 3.25021005, 3.30490994, 3.35956001, 3.41415000,
    3.46868992, 3.52316999, 3.57758999, 3.63195992, 3.68628001, 3.74054003, 3.79475999, 3.84892988,
    3.90302992, 3.95708990, 4.01110983, 4.06506014, 4.11896992,
    3.08008003, 3.13491988, 3.18969989, 3.24440002, 3.29905009, 3.35363007, 3.40814996, 3.46263003,
    3.51704001, 3.57140994, 3.62572002, 3.67998004, 3.73417997, 3.78834009, 3.84243989, 3.89650011,
    3.95051003, 4.00446987, 4.05837011, 4.11223984, 4.16604996,
    3.12948990, 3.18422008, 3.23887992, 3.29346991, 3.34800005, 3.40247011, 3.45688009, 3.51124001,
    3.56554008, 3.61980009, 3.67400002, 3.72814989, 3.78224993, 3.83628988, 3.89030004, 3.94425011,
    3.99815989, 4.05203009, 4.10583019, 4.15960979, 4.21332979,
    3.17898989, 3.23360991, 3.28815007, 3.34263992, 3.39705992, 3.45142007, 3.50570989, 3.55996990,
    3.61416006, 3.66830993, 3.72240996, 3.77644992, 3.83046007, 3.88439989, 3.93830991, 3.99216008,
    4.04597998, 4.09974003, 4.15347004, 4.20714998, 4.26077986,
    3.22854996, 3.28307009, 3.33751011, 3.39188004, 3.44619989, 3.50045991, 3.55466008, 3.60879993,
    3.66289997, 3.71693993, 3.77095008, 3.82488990, 3.87879992, 3.93265009, 3.98645997, 4.04021978,
    4.09394979, 4.14762020, 4.20126009, 4.25484991, 4.30840015 };
  neck_max_separation.putHost(sep_buffer);
  sp_neck_max_separation.putHost(std::vector<float>(sep_buffer.begin(), sep_buffer.end()));

  // Load the values table for neck GB
  std::vector<double> val_buffer = {
    0.03815110, 0.03385870, 0.03017760, 0.02700300, 0.02425060, 0.02185290, 0.01975470, 0.01791090,
    0.01628440, 0.01484420, 0.01356470, 0.01242430, 0.01140470, 0.01049060, 0.00966876, 0.00892800,
    0.00825870, 0.00765255, 0.00710237, 0.00660196, 0.00614589,
    0.03961980, 0.03518370, 0.03137670, 0.02809110, 0.02524090, 0.02275630, 0.02058080, 0.01866810,
    0.01697990, 0.01548430, 0.01415500, 0.01296960, 0.01190940, 0.01095840, 0.01010310, 0.00933189,
    0.00863480, 0.00800326, 0.00742986, 0.00690814, 0.00643255,
    0.04104800, 0.03647380, 0.03254560, 0.02915320, 0.02620840, 0.02363990, 0.02138970, 0.01941020,
    0.01766220, 0.01611290, 0.01473510, 0.01350590, 0.01240610, 0.01141920, 0.01053120, 0.00973027,
    0.00900602, 0.00834965, 0.00775350, 0.00721091, 0.00671609,
    0.04243650, 0.03772950, 0.03368460, 0.03018930, 0.02715330, 0.02450380, 0.02218130, 0.02013710,
    0.01833100, 0.01672950, 0.01530470, 0.01403300, 0.01289460, 0.01187270, 0.01095290, 0.01012290,
    0.00937212, 0.00869147, 0.00807306, 0.00751003, 0.00699641,
    0.04378610, 0.03895160, 0.03479440, 0.03119980, 0.02807580, 0.02534790, 0.02295550, 0.02084870,
    0.01898640, 0.01733430, 0.01586370, 0.01455070, 0.01337480, 0.01231880, 0.01136790, 0.01050960,
    0.00973290, 0.00902853, 0.00838835, 0.00780533, 0.00727330,
    0.04509790, 0.04014060, 0.03587530, 0.03218510, 0.02897610, 0.02617260, 0.02371250, 0.02154510,
    0.01962820, 0.01792700, 0.01641210, 0.01505880, 0.01384650, 0.01275730, 0.01177610, 0.01089020,
    0.01008820, 0.00936068, 0.00869923, 0.00809665, 0.00754661,
    0.04637290, 0.04129760, 0.03692810, 0.03314560, 0.02985470, 0.02697800, 0.02445250, 0.02222640,
    0.02025670, 0.01850780, 0.01694980, 0.01555750, 0.01430960, 0.01318810, 0.01217750, 0.01126460,
    0.01043800, 0.00968781, 0.00900559, 0.00838388, 0.00781622,
    0.04761230, 0.04242330, 0.03795340, 0.03408200, 0.03071180, 0.02776450, 0.02517570, 0.02289270,
    0.02087180, 0.01907670, 0.01747680, 0.01604660, 0.01476420, 0.01361120, 0.01257190, 0.01163280,
    0.01078210, 0.01000990, 0.00930735, 0.00866695, 0.00808206,
    0.04881710, 0.04351860, 0.03895200, 0.03499470, 0.03154810, 0.02853240, 0.02588240, 0.02354430,
    0.02147380, 0.01963390, 0.01799340, 0.01652620, 0.01521030, 0.01402670, 0.01295950, 0.01199470,
    0.01112060, 0.01032680, 0.00960445, 0.00894579, 0.00834405,
    0.04998830, 0.04458450, 0.03992460, 0.03588440, 0.03236400, 0.02928220, 0.02657290, 0.02418150,
    0.02206290, 0.02017940, 0.01849940, 0.01699640, 0.01564790, 0.01443450, 0.01334010, 0.01235040,
    0.01145340, 0.01063860, 0.00989687, 0.00922037, 0.00860216,
    0.05112720, 0.04562190, 0.04087200, 0.03675180, 0.03315990, 0.03001420, 0.02724750, 0.02480450,
    0.02263920, 0.02071350, 0.01899520, 0.01745740, 0.01607710, 0.01483480, 0.01371380, 0.01269980,
    0.01178050, 0.01094520, 0.01018460, 0.00949067, 0.00885636,
    0.05223480, 0.04663150, 0.04179480, 0.03759730, 0.03393650, 0.03072900, 0.02790670, 0.02541360,
    0.02320300, 0.02123630, 0.01948090, 0.01790920, 0.01649800, 0.01522750, 0.01408070, 0.01304300,
    0.01210200, 0.01124660, 0.01046760, 0.00975668, 0.00910664,
    0.05331230, 0.04761450, 0.04269400, 0.03842180, 0.03469420, 0.03142680, 0.02855070, 0.02600900,
    0.02375470, 0.02174820, 0.01995660, 0.01835200, 0.01691080, 0.01561280, 0.01444080, 0.01338010,
    0.01241790, 0.01154300, 0.01074600, 0.01001840, 0.00935302,
    0.05436060, 0.04857160, 0.04357000, 0.03922570, 0.03543350, 0.03210820, 0.02918000, 0.02659130,
    0.02429430, 0.02224920, 0.02042250, 0.01878590, 0.01731550, 0.01599080, 0.01479430, 0.01371110,
    0.01272820, 0.01183430, 0.01101970, 0.01027590, 0.00959549,
    0.05538070, 0.04950370, 0.04442390, 0.04000970, 0.03615510, 0.03277360, 0.02979490, 0.02716050,
    0.02482220, 0.02273960, 0.02087880, 0.01921110, 0.01771220, 0.01636150, 0.01514130, 0.01403610,
    0.01303300, 0.01212060, 0.01128880, 0.01052920, 0.00983409,
    0.05637380, 0.05041160, 0.04525620, 0.04077450, 0.03685930, 0.03342350, 0.03039580, 0.02771710,
    0.02533870, 0.02321970, 0.02132570, 0.01962770, 0.01810130, 0.01672520, 0.01548170, 0.01435520,
    0.01333250, 0.01240190, 0.01155340, 0.01077830, 0.01006880,
    0.05734060, 0.05129630, 0.04606760, 0.04152060, 0.03754680, 0.03405830, 0.03098300, 0.02826140,
    0.02584410, 0.02368960, 0.02176340, 0.02003600, 0.01848260, 0.01708200, 0.01581580, 0.01466850,
    0.01362660, 0.01267830, 0.01181350, 0.01102320, 0.01029980,
    0.05828220, 0.05215840, 0.04685890, 0.04224860, 0.03821800, 0.03467840, 0.03155710, 0.02879380,
    0.02633860, 0.02414970, 0.02219220, 0.02043620, 0.01885660, 0.01743190, 0.01614370, 0.01497610,
    0.01391540, 0.01294990, 0.01206910, 0.01126410, 0.01052690,
    0.05919940, 0.05299870, 0.04763070, 0.04295900, 0.03887340, 0.03528430, 0.03211820, 0.02931440,
    0.02682250, 0.02460020, 0.02261210, 0.02082830, 0.01922320, 0.01777510, 0.01646540, 0.01527800,
    0.01419910, 0.01321670, 0.01232040, 0.01150090, 0.01075040,
    0.06009320, 0.05381800, 0.04838360, 0.04365250, 0.03951360, 0.03587640, 0.03266690, 0.02982370,
    0.02729610, 0.02504130, 0.02302360, 0.02121260, 0.01958260, 0.01811180, 0.01678110, 0.01557440,
    0.01447780, 0.01347890, 0.01256730, 0.01173380, 0.01097020,
    0.06096420, 0.05461690, 0.04911830, 0.04432950, 0.04013880, 0.03645500, 0.03320330, 0.03032200,
    0.02775960, 0.02547320, 0.02342660, 0.02158920, 0.01993510, 0.01844200, 0.01709090, 0.01586540,
    0.01475140, 0.01373650, 0.01281010, 0.01196270, 0.01118630 };
  neck_max_value.putHost(val_buffer);
  sp_neck_max_value.putHost(std::vector<float>(val_buffer.begin(), val_buffer.end()));

  // Upload data to the device, if applicable
#ifdef STORMM_USE_HPC
  neck_max_separation.upload();
  neck_max_value.upload();
  sp_neck_max_separation.upload();
  sp_neck_max_value.upload();
#endif
}

//-------------------------------------------------------------------------------------------------
NeckGeneralizedBornTable::NeckGeneralizedBornTable(const int table_size_in,
                                                   const double neck_cut_in,
                                                   const double kscale_in,
                                                   const std::vector<double> &max_sep_in,
                                                   const std::vector<double> &max_val_in) :
    table_size{table_size_in}, neck_cut{neck_cut_in}, kscale{kscale_in},
    neck_max_separation{static_cast<size_t>(table_size * table_size), "gb_table_maxsep"},
    neck_max_value{static_cast<size_t>(table_size * table_size), "gb_table_maxval"},
    sp_neck_max_separation{static_cast<size_t>(table_size * table_size), "gb_table_maxsep"},
    sp_neck_max_value{static_cast<size_t>(table_size * table_size), "gb_table_maxval"}
{
  if (static_cast<int>(max_sep_in.size()) != table_size_in * table_size_in) {
    rtErr("\"Neck\" Generalized Born tables cannot be assembled with size " +
          std::to_string(table_size_in) + " with data of length " +
          std::to_string(max_sep_in.size()) + ".", "NeckGeneralizedBornTable");
  }
  if (max_sep_in.size() != max_val_in.size()) {
    rtErr("Tables cannot be of unequal sizes (" + std::to_string(max_sep_in.size()) + " and " +
          std::to_string(max_val_in.size()) + ").", "NeckGeneralizedBornTable");
  }
  neck_max_separation.putHost(max_sep_in);
  neck_max_value.putHost(max_val_in);
  sp_neck_max_separation.putHost(std::vector<float>(max_sep_in.begin(), max_sep_in.end()));
  sp_neck_max_value.putHost(std::vector<float>(max_val_in.begin(), max_val_in.end()));

  // Upload data to the device, if applicable
#ifdef STORMM_USE_HPC
  neck_max_separation.upload();
  neck_max_value.upload();
  sp_neck_max_separation.upload();
  sp_neck_max_value.upload();
#endif
}

//-------------------------------------------------------------------------------------------------
int NeckGeneralizedBornTable::size() const {
  return table_size;
}

//-------------------------------------------------------------------------------------------------
double NeckGeneralizedBornTable::getMaxSeparation(int i_type, int j_type) const {
  return neck_max_separation.readHost((j_type * table_size) + i_type);
}

//-------------------------------------------------------------------------------------------------
double NeckGeneralizedBornTable::getMaxValue(int i_type, int j_type) const {
  return neck_max_value.readHost((j_type * table_size) + i_type);
}

//-------------------------------------------------------------------------------------------------
const NeckGeneralizedBornKit<double>
NeckGeneralizedBornTable::dpData(const HybridTargetLevel tier) const {
  return NeckGeneralizedBornKit<double>(table_size, neck_cut, kscale,
                                        neck_max_separation.data(tier),
                                        neck_max_value.data(tier));
}

//-------------------------------------------------------------------------------------------------
const NeckGeneralizedBornKit<float>
NeckGeneralizedBornTable::spData(const HybridTargetLevel tier) const {
  return NeckGeneralizedBornKit<float>(table_size, neck_cut, kscale,
                                       sp_neck_max_separation.data(tier),
                                       sp_neck_max_value.data(tier));
}

//-------------------------------------------------------------------------------------------------
const NeckGeneralizedBornTable* NeckGeneralizedBornTable::getSelfPointer() const {
  return this;
}

} // namespace generalized_born_defaults
} // namespace stormm

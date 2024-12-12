//!
//!  Everything related to the $\mathrm{RotSum}$ operation.
//!

use ark_ff::{Field, Zero};
use lattirust_ring::{Cyclotomic, PolyRing};

use crate::ark_base::*;
use crate::rings::SuitableRing;

/// An implementation of the $\mathrm{RotSum}$ operation from lemma 2.1 of the Latticefold paper.
///
/// The formula is
/// $$
/// \mathrm{RotSum}(a\in \mathcal{R}\_p, \vec{B}\in\mathbb{Z}\_{p^\tau}^d)=\sum\_{i=1}^d B_i\cdot
///     \mathrm{Coeff}(X^{i-1}a)\in\mathbb{Z}\_{p^\tau}^d.
/// $$
///
/// The function computes the linear combination of successive rotations of a ring element
/// (treated as their coefficient vectors) with coefficients from the base field of the NTT form
/// of the ring.
///
/// The function iterates through each element `b_i` of the input slice `b`
/// (of length `R::CoefficientRepresentation::dimension()`), computes the
/// corresponding rotated coefficients from `a`, and multiplies each coefficient
/// by the corresponding `b_i` from the base ring. These products are then summed
/// together across all coefficients for each element of the base ring.
///
/// This function requires that the length of `b` matches the dimension of the
/// `R::CoefficientRepresentation`, as enforced by the `assert_eq!` check.
///
/// # Parameters:
/// - `a`: A coefficient representation of a ring element.
/// - `b`: A slice of elements from the base ring `R::BaseRing`
///   (note that for a suitable ring it is guaranteed to be a field) that will be used
///   to scale the coefficients from `a`.
///
/// # Returns:
/// A `Vec<R::BaseRing>` representing the resulting sum after multiplying and
/// summing the scaled coefficients.
///
/// # Panics:
/// The function panics if the length of `b` does not match the degree of
/// of the cyclotomic polynomial (which in this cases is equal to `R::CoefficientRepresentation::dimension()`).
///
pub fn rot_sum<R: SuitableRing>(
    a: R::CoefficientRepresentation,
    b: &[R::BaseRing],
) -> Vec<R::BaseRing> {
    assert_eq!(b.len(), R::CoefficientRepresentation::dimension());

    let mut acc = vec![R::BaseRing::zero(); R::CoefficientRepresentation::dimension()];

    for (b_i, x_i_a) in b.iter().zip(a.into_rot_iter()) {
        for (acc_j, x) in acc.iter_mut().zip(x_i_a.into_coeffs()) {
            *acc_j += <R::BaseRing as Field>::from_base_prime_field(x) * b_i;
        }
    }

    acc
}

/// Computes the sum $\sum\_{i=1}^{2k} \mathrm{RotSum}(\rho\_i, \mathrm{NTT}(\theta\_i))$.
///
/// The function takes two slices: one of coefficient representations of the rho's (`rho_s`),
/// and one of the theta vectors (`theta_s`) from the base field. For each pair of elements
/// `(rho_i, theta_i)` from `rho_s` and `theta_s`, it computes the rotated sums
/// of coefficients from `rho_s` and vectors from `theta_i`. The results of these
/// rotated sums are then summed together.
///
/// # Parameters:
/// - `rho_s`: A slice of ring elements in the coefficient representation `R::CoefficientRepresentation`.
///   The length of `rho_s` must match the length of `theta_s`.
/// - `theta_s`: A slice of vectors of ring elements in the NTT form (`Vec<R>`).
///
/// # Returns:
/// A `Vec<R>` which is the sum of the rotated sums $\mathrm{RotSum}(\rho\_i, \mathrm{NTT}(\theta\_i))$.
///
/// # Panics:
/// The function will panic if the lengths of `rho_s` and `theta_s` do not match,
/// as enforced by the `assert_eq!` check. It will also panic if the final result
/// cannot be promoted from coefficients due to dimensional or capacity issues,
/// as indicated by the `expect` call in `R::promote_from_coeffs`.
///
pub fn rot_lin_combination<R: SuitableRing>(
    rho_s: &[R::CoefficientRepresentation],
    theta_s: &[Vec<R>],
) -> Vec<R> {
    assert_eq!(rho_s.len(), theta_s.len());

    // Here we assume that `R::flatten_to_coeffs` transforms the capacity of a vector as `cap * R::dimension()`.
    let mut res = R::flatten_to_coeffs(vec![R::zero(); theta_s[0].len()]);

    for (&rho_i, theta_i) in rho_s.iter().zip(theta_s) {
        let theta_flat = R::flatten_to_coeffs(theta_i.clone());

        let sum = rot_sum::<R>(rho_i, &theta_flat);

        for (i, x) in sum.iter().enumerate() {
            res[i] += x;
        }
    }

    R::promote_from_coeffs(res).expect("The length and the capacity should be correct")
}

#[cfg(test)]
mod tests {
    use ark_ff::UniformRand;
    use lattirust_ring::cyclotomic_ring::models::goldilocks::{Fq, Fq3};

    use super::*;
    use crate::rings::{GoldilocksRingNTT, GoldilocksRingPoly};

    #[test]
    fn test_rot_sum_with_coeffs() {
        let mut rng = ark_std::test_rng();
        let a = GoldilocksRingPoly::rand(&mut rng);
        let b = GoldilocksRingPoly::rand(&mut rng);

        // RotSum(a, coeff(b)) = coeff(a * b)
        assert_eq!(
            rot_sum::<GoldilocksRingNTT>(
                a,
                &b.coeffs()
                    .iter()
                    .map(|&x| Fq3::from_base_prime_field(x))
                    .collect::<Vec<Fq3>>()
            ),
            (a * b)
                .into_coeffs()
                .into_iter()
                .map(Fq3::from_base_prime_field)
                .collect::<Vec<Fq3>>()
        );
    }

    #[test]
    fn test_rot_sum_with_ring_elems() {
        let mut rng = ark_std::test_rng();

        let a = GoldilocksRingPoly::rand(&mut rng);
        let b: Vec<GoldilocksRingPoly> = (0..Fq3::extension_degree())
            .map(|_| GoldilocksRingPoly::rand(&mut rng))
            .collect();

        let capital_b: Vec<Fq3> = (0..GoldilocksRingPoly::dimension())
            .map(|i| {
                Fq3::from_base_prime_field_elems([
                    b[0].coeffs()[i],
                    b[1].coeffs()[i],
                    b[2].coeffs()[i],
                ])
                .unwrap()
            })
            .collect();

        // RotSum(a, B) = \sum coeff(a * b_i) * Y^{i-1}
        assert_eq!(
            rot_sum::<GoldilocksRingNTT>(a, &capital_b),
            (0..GoldilocksRingPoly::dimension())
                .map(|i| {
                    Fq3::from_base_prime_field_elems([
                        (a * b[0]).coeffs()[i],
                        (a * b[1]).coeffs()[i],
                        (a * b[2]).coeffs()[i],
                    ])
                    .unwrap()
                })
                .collect::<Vec<Fq3>>()
        );
    }

    #[test]
    fn test_rot_lin_combination() {
        let rho_s = vec![
            GoldilocksRingPoly::from(vec![
                Fq::from(3176501537275927775u64),
                Fq::from(4431515930232641055u64),
                Fq::from(13221694751451514597u64),
                Fq::from(5465102077112450573u64),
                Fq::from(7182954176643397363u64),
                Fq::from(16712294236572839616u64),
                Fq::from(2067490596382165132u64),
                Fq::from(10217230022534865163u64),
                Fq::from(6229993698204923738u64),
                Fq::from(1325193567257335506u64),
                Fq::from(13078047754975611428u64),
                Fq::from(4834528168480107204u64),
                Fq::from(17612599349192843794u64),
                Fq::from(15747147843029123525u64),
                Fq::from(1554538465021420955u64),
                Fq::from(4600541983470645889u64),
                Fq::from(14439397527338063444u64),
                Fq::from(16713363142778527590u64),
                Fq::from(11812384261521111838u64),
                Fq::from(17448028986017570869u64),
                Fq::from(10865383156233274881u64),
                Fq::from(3383399182301445293u64),
                Fq::from(10936029056371001221u64),
                Fq::from(15242877653219303325u64),
            ]),
            GoldilocksRingPoly::from(vec![
                Fq::from(12772015572104235724u64),
                Fq::from(3038266509051673063u64),
                Fq::from(830468003174106144u64),
                Fq::from(1657201770232657737u64),
                Fq::from(17632306365478334080u64),
                Fq::from(16532585872688629458u64),
                Fq::from(3901882944239384404u64),
                Fq::from(17631340313161239038u64),
                Fq::from(608476189009328641u64),
                Fq::from(11056657407656114942u64),
                Fq::from(11102497970140551326u64),
                Fq::from(12464256964897601353u64),
                Fq::from(269748484779416232u64),
                Fq::from(10243040561263411480u64),
                Fq::from(13668066756487173204u64),
                Fq::from(10915381335733689310u64),
                Fq::from(4244243121477916167u64),
                Fq::from(3832142791392918443u64),
                Fq::from(316851653036921491u64),
                Fq::from(16234292544561705006u64),
                Fq::from(2029597689787369711u64),
                Fq::from(5388502217815959971u64),
                Fq::from(1914802618499239668u64),
                Fq::from(12153153672757460902u64),
            ]),
            GoldilocksRingPoly::from(vec![
                Fq::from(17393762974846229185u64),
                Fq::from(16744594376952940276u64),
                Fq::from(18245951515647764122u64),
                Fq::from(3685771524837185649u64),
                Fq::from(277274257190934084u64),
                Fq::from(9835391571776287244u64),
                Fq::from(2069820602651091944u64),
                Fq::from(8450389677968116299u64),
                Fq::from(4054156578124104751u64),
                Fq::from(18211771542183737663u64),
                Fq::from(7142882814424378912u64),
                Fq::from(1742921454545732202u64),
                Fq::from(15547846078353729808u64),
                Fq::from(1550252619098832286u64),
                Fq::from(3871471667468774288u64),
                Fq::from(9178691481501245860u64),
                Fq::from(1014750535827442u64),
                Fq::from(12954068300922993304u64),
                Fq::from(11652068639258192282u64),
                Fq::from(104271770142799833u64),
                Fq::from(3352524758668218123u64),
                Fq::from(12800226862690544761u64),
                Fq::from(18423001377572820480u64),
                Fq::from(6594862545730709887u64),
            ]),
        ];

        let theta_s = vec![
            vec![
                GoldilocksRingNTT::from(vec![
                    Fq3::new(
                        Fq::from(10657762469820564374u64),
                        Fq::from(13124963420535130615u64),
                        Fq::from(4009228708511043147u64),
                    ),
                    Fq3::new(
                        Fq::from(14305463036021180822u64),
                        Fq::from(12496786129003557812u64),
                        Fq::from(11589867157344105204u64),
                    ),
                    Fq3::new(
                        Fq::from(15389545485034963067u64),
                        Fq::from(8264739475908697585u64),
                        Fq::from(8992666529373462528u64),
                    ),
                    Fq3::new(
                        Fq::from(5255545208733385525u64),
                        Fq::from(12203460203780646213u64),
                        Fq::from(866129832992083504u64),
                    ),
                    Fq3::new(
                        Fq::from(6537292673425061051u64),
                        Fq::from(298327783385136554u64),
                        Fq::from(17044980290698036648u64),
                    ),
                    Fq3::new(
                        Fq::from(13910712063845964384u64),
                        Fq::from(15508196130543392522u64),
                        Fq::from(16628506396074271880u64),
                    ),
                    Fq3::new(
                        Fq::from(7949006387678825284u64),
                        Fq::from(1104655596793289123u64),
                        Fq::from(7512076554559487881u64),
                    ),
                    Fq3::new(
                        Fq::from(13244979870851355903u64),
                        Fq::from(14568618605873407212u64),
                        Fq::from(6572380305465045567u64),
                    ),
                ]),
                GoldilocksRingNTT::from(vec![
                    Fq3::new(
                        Fq::from(11870588123800510954u64),
                        Fq::from(9946741634701123808u64),
                        Fq::from(4844329685037136431u64),
                    ),
                    Fq3::new(
                        Fq::from(15341191814094103905u64),
                        Fq::from(1552677365453658370u64),
                        Fq::from(17122433861781042570u64),
                    ),
                    Fq3::new(
                        Fq::from(17759918033711952814u64),
                        Fq::from(2743353454026887314u64),
                        Fq::from(8636188326248148928u64),
                    ),
                    Fq3::new(
                        Fq::from(15514752800134252955u64),
                        Fq::from(3809218950601361321u64),
                        Fq::from(12312963399730956268u64),
                    ),
                    Fq3::new(
                        Fq::from(5838488256440040115u64),
                        Fq::from(170586405955424791u64),
                        Fq::from(16225200886416819310u64),
                    ),
                    Fq3::new(
                        Fq::from(4516619279923596297u64),
                        Fq::from(10625667863499443133u64),
                        Fq::from(6806052842469481847u64),
                    ),
                    Fq3::new(
                        Fq::from(16495583648158152584u64),
                        Fq::from(5133514671042811519u64),
                        Fq::from(13501107762810765614u64),
                    ),
                    Fq3::new(
                        Fq::from(8556538166913825858u64),
                        Fq::from(13645963490328782150u64),
                        Fq::from(16185848549638148805u64),
                    ),
                ]),
                GoldilocksRingNTT::from(vec![
                    Fq3::new(
                        Fq::from(15242515788959265127u64),
                        Fq::from(11484430514333675184u64),
                        Fq::from(4645629102841269342u64),
                    ),
                    Fq3::new(
                        Fq::from(1420906651893672818u64),
                        Fq::from(14511527686584485704u64),
                        Fq::from(4019368467235939780u64),
                    ),
                    Fq3::new(
                        Fq::from(6432841109420381910u64),
                        Fq::from(2379970813695579460u64),
                        Fq::from(16108004456390486413u64),
                    ),
                    Fq3::new(
                        Fq::from(13399550195234911174u64),
                        Fq::from(6388020423213188602u64),
                        Fq::from(9594868308224470184u64),
                    ),
                    Fq3::new(
                        Fq::from(7920868308701161530u64),
                        Fq::from(5276946360535732918u64),
                        Fq::from(17574966081460745192u64),
                    ),
                    Fq3::new(
                        Fq::from(11342097745811234830u64),
                        Fq::from(14425103415146572053u64),
                        Fq::from(503224155324747872u64),
                    ),
                    Fq3::new(
                        Fq::from(2808235238362036791u64),
                        Fq::from(5005263963108873187u64),
                        Fq::from(16613197613331748953u64),
                    ),
                    Fq3::new(
                        Fq::from(11902766678878194632u64),
                        Fq::from(4183991245447786665u64),
                        Fq::from(1618047714598681981u64),
                    ),
                ]),
            ],
            vec![
                GoldilocksRingNTT::from(vec![
                    Fq3::new(
                        Fq::from(12115513894979888991u64),
                        Fq::from(831492639124919797u64),
                        Fq::from(5099924383830363545u64),
                    ),
                    Fq3::new(
                        Fq::from(8334164740025221756u64),
                        Fq::from(14277555808239109159u64),
                        Fq::from(14414977308918277632u64),
                    ),
                    Fq3::new(
                        Fq::from(2802144098390705189u64),
                        Fq::from(15560520802099246211u64),
                        Fq::from(9622511342221168212u64),
                    ),
                    Fq3::new(
                        Fq::from(16368574357912424107u64),
                        Fq::from(2860613981976064479u64),
                        Fq::from(4771949792573268947u64),
                    ),
                    Fq3::new(
                        Fq::from(5244795760798913513u64),
                        Fq::from(8123967081328510661u64),
                        Fq::from(14578430994012767928u64),
                    ),
                    Fq3::new(
                        Fq::from(5033077145447971279u64),
                        Fq::from(546532940620769021u64),
                        Fq::from(16213800138522747930u64),
                    ),
                    Fq3::new(
                        Fq::from(16518446145116382763u64),
                        Fq::from(11893632984880492997u64),
                        Fq::from(3760330906359743252u64),
                    ),
                    Fq3::new(
                        Fq::from(10201857189930300557u64),
                        Fq::from(7470120683054760655u64),
                        Fq::from(10527800928261867971u64),
                    ),
                ]),
                GoldilocksRingNTT::from(vec![
                    Fq3::new(
                        Fq::from(11370241597087230387u64),
                        Fq::from(10290746092337449848u64),
                        Fq::from(414303224995321709u64),
                    ),
                    Fq3::new(
                        Fq::from(16036965224018329056u64),
                        Fq::from(10616732916918694076u64),
                        Fq::from(13530953995307178796u64),
                    ),
                    Fq3::new(
                        Fq::from(3318856251365656862u64),
                        Fq::from(6934703379937920412u64),
                        Fq::from(17633083368184917936u64),
                    ),
                    Fq3::new(
                        Fq::from(4673716009888813010u64),
                        Fq::from(3153404660295176257u64),
                        Fq::from(13237192770515254077u64),
                    ),
                    Fq3::new(
                        Fq::from(13666750031925975117u64),
                        Fq::from(3414857594022568181u64),
                        Fq::from(7762116978150506244u64),
                    ),
                    Fq3::new(
                        Fq::from(6941716297597842418u64),
                        Fq::from(201347344700946031u64),
                        Fq::from(8796306432393489059u64),
                    ),
                    Fq3::new(
                        Fq::from(2826121039681926212u64),
                        Fq::from(14406062495830747443u64),
                        Fq::from(6787050749566664671u64),
                    ),
                    Fq3::new(
                        Fq::from(4267934789887714727u64),
                        Fq::from(8184440252650452322u64),
                        Fq::from(7063331600913089197u64),
                    ),
                ]),
                GoldilocksRingNTT::from(vec![
                    Fq3::new(
                        Fq::from(1057619362498405958u64),
                        Fq::from(14132593295000936311u64),
                        Fq::from(778288501325222844u64),
                    ),
                    Fq3::new(
                        Fq::from(5156252787319925221u64),
                        Fq::from(2005631339388013655u64),
                        Fq::from(6800683298592284644u64),
                    ),
                    Fq3::new(
                        Fq::from(12435812410644450340u64),
                        Fq::from(13807734065384196529u64),
                        Fq::from(12618804004437426850u64),
                    ),
                    Fq3::new(
                        Fq::from(7450108733651691952u64),
                        Fq::from(12595045956724047917u64),
                        Fq::from(1228204942079097619u64),
                    ),
                    Fq3::new(
                        Fq::from(16016377034106130664u64),
                        Fq::from(9738927012705256472u64),
                        Fq::from(15800628473620236772u64),
                    ),
                    Fq3::new(
                        Fq::from(13144933341810616865u64),
                        Fq::from(12599599890183199476u64),
                        Fq::from(11629560998570673175u64),
                    ),
                    Fq3::new(
                        Fq::from(16462442903005482401u64),
                        Fq::from(4881714654374311256u64),
                        Fq::from(2763207197959151835u64),
                    ),
                    Fq3::new(
                        Fq::from(5236783854299915347u64),
                        Fq::from(8290795181217338883u64),
                        Fq::from(1848465477998056285u64),
                    ),
                ]),
            ],
            vec![
                GoldilocksRingNTT::from(vec![
                    Fq3::new(
                        Fq::from(12140311057295504440u64),
                        Fq::from(10595884708328810007u64),
                        Fq::from(13593226595043988994u64),
                    ),
                    Fq3::new(
                        Fq::from(13436656350607313037u64),
                        Fq::from(12229892596877328739u64),
                        Fq::from(9319801939418154159u64),
                    ),
                    Fq3::new(
                        Fq::from(12420817867601258274u64),
                        Fq::from(12666487015478164843u64),
                        Fq::from(2193938772994168660u64),
                    ),
                    Fq3::new(
                        Fq::from(15849450528139170991u64),
                        Fq::from(13832774759315991642u64),
                        Fq::from(3408734832133577312u64),
                    ),
                    Fq3::new(
                        Fq::from(13099559230793721491u64),
                        Fq::from(12901171124983892608u64),
                        Fq::from(17248733471520632964u64),
                    ),
                    Fq3::new(
                        Fq::from(12719575413580328783u64),
                        Fq::from(12842983474236239556u64),
                        Fq::from(16931631236564036085u64),
                    ),
                    Fq3::new(
                        Fq::from(5645831881453102036u64),
                        Fq::from(5342521708667963235u64),
                        Fq::from(7714619877603915774u64),
                    ),
                    Fq3::new(
                        Fq::from(1298699196556602978u64),
                        Fq::from(6432856332967067723u64),
                        Fq::from(12103290744719120900u64),
                    ),
                ]),
                GoldilocksRingNTT::from(vec![
                    Fq3::new(
                        Fq::from(14192596637610880286u64),
                        Fq::from(302603983284521212u64),
                        Fq::from(1554091430383766001u64),
                    ),
                    Fq3::new(
                        Fq::from(15235589363489599773u64),
                        Fq::from(6361068694970402894u64),
                        Fq::from(11351496496142666695u64),
                    ),
                    Fq3::new(
                        Fq::from(17798936219494665236u64),
                        Fq::from(17962457795929401978u64),
                        Fq::from(2393501487684627430u64),
                    ),
                    Fq3::new(
                        Fq::from(11633635576539857675u64),
                        Fq::from(6820539544375496400u64),
                        Fq::from(7026687734296337460u64),
                    ),
                    Fq3::new(
                        Fq::from(3264513787526822383u64),
                        Fq::from(7513232642957019192u64),
                        Fq::from(10197468329800862154u64),
                    ),
                    Fq3::new(
                        Fq::from(9295087866276072263u64),
                        Fq::from(3923753298006276885u64),
                        Fq::from(17188588220609937953u64),
                    ),
                    Fq3::new(
                        Fq::from(15618300783890311120u64),
                        Fq::from(18282329808415228385u64),
                        Fq::from(12853250083733216116u64),
                    ),
                    Fq3::new(
                        Fq::from(2535369068885733448u64),
                        Fq::from(2951896220056883455u64),
                        Fq::from(14198387874427779401u64),
                    ),
                ]),
                GoldilocksRingNTT::from(vec![
                    Fq3::new(
                        Fq::from(9252710647055088837u64),
                        Fq::from(16619479957846686549u64),
                        Fq::from(7078168544294434064u64),
                    ),
                    Fq3::new(
                        Fq::from(285579269951443452u64),
                        Fq::from(12636207664827152311u64),
                        Fq::from(4837517500951888315u64),
                    ),
                    Fq3::new(
                        Fq::from(3426615841556639802u64),
                        Fq::from(3695048890694829804u64),
                        Fq::from(3066505982824881996u64),
                    ),
                    Fq3::new(
                        Fq::from(15138353857325256724u64),
                        Fq::from(16029794379726145376u64),
                        Fq::from(12218332594708140109u64),
                    ),
                    Fq3::new(
                        Fq::from(8561470032224929624u64),
                        Fq::from(17405176291630133625u64),
                        Fq::from(5745530902628363247u64),
                    ),
                    Fq3::new(
                        Fq::from(10692454097795272437u64),
                        Fq::from(6654090244064764695u64),
                        Fq::from(6886701096438338327u64),
                    ),
                    Fq3::new(
                        Fq::from(5279718514473932197u64),
                        Fq::from(3903201918770756414u64),
                        Fq::from(4055677395422880309u64),
                    ),
                    Fq3::new(
                        Fq::from(10483642917291729686u64),
                        Fq::from(6483375199820341235u64),
                        Fq::from(5032927012063871650u64),
                    ),
                ]),
            ],
        ];

        let res = rot_lin_combination(&rho_s, &theta_s);

        let expected = vec![
            GoldilocksRingNTT::from(vec![
                Fq3::new(
                    Fq::from(14333051603278409827u64),
                    Fq::from(3452210451626351104u64),
                    Fq::from(11553319917643724446u64),
                ),
                Fq3::new(
                    Fq::from(4271198654290429726u64),
                    Fq::from(5457779707607808166u64),
                    Fq::from(4492361609639224251u64),
                ),
                Fq3::new(
                    Fq::from(18403103392528269801u64),
                    Fq::from(10748780410942685903u64),
                    Fq::from(4825826961141717824u64),
                ),
                Fq3::new(
                    Fq::from(6968286132664619352u64),
                    Fq::from(8248545517236422401u64),
                    Fq::from(14619073437634613110u64),
                ),
                Fq3::new(
                    Fq::from(6494468683726654237u64),
                    Fq::from(6736819944922364764u64),
                    Fq::from(7933269282055173624u64),
                ),
                Fq3::new(
                    Fq::from(1411312866761811864u64),
                    Fq::from(11740388575380027346u64),
                    Fq::from(12440697690963923599u64),
                ),
                Fq3::new(
                    Fq::from(11430892440023762418u64),
                    Fq::from(8733019049023381480u64),
                    Fq::from(12142894848969621232u64),
                ),
                Fq3::new(
                    Fq::from(5304908433760179753u64),
                    Fq::from(2683457664236894477u64),
                    Fq::from(4277151749528966082u64),
                ),
            ]),
            GoldilocksRingNTT::from(vec![
                Fq3::new(
                    Fq::from(12865678513040868880u64),
                    Fq::from(8330946251785883254u64),
                    Fq::from(13727156933229298000u64),
                ),
                Fq3::new(
                    Fq::from(5764046672888055774u64),
                    Fq::from(5899144915470388003u64),
                    Fq::from(11217054836791419852u64),
                ),
                Fq3::new(
                    Fq::from(541021512657825698u64),
                    Fq::from(14183758266840145326u64),
                    Fq::from(10360899138354881232u64),
                ),
                Fq3::new(
                    Fq::from(11328191143048337580u64),
                    Fq::from(15950133360099954580u64),
                    Fq::from(2565741296849987424u64),
                ),
                Fq3::new(
                    Fq::from(12497879903222337532u64),
                    Fq::from(7978382415323963175u64),
                    Fq::from(13389595069071535811u64),
                ),
                Fq3::new(
                    Fq::from(2555728133273857737u64),
                    Fq::from(12506608867240787524u64),
                    Fq::from(201605385374345405u64),
                ),
                Fq3::new(
                    Fq::from(16412292333012688894u64),
                    Fq::from(12828121542832933949u64),
                    Fq::from(9253892227253570463u64),
                ),
                Fq3::new(
                    Fq::from(16722100444222153203u64),
                    Fq::from(18107828994824775518u64),
                    Fq::from(9546850266634560739u64),
                ),
            ]),
            GoldilocksRingNTT::from(vec![
                Fq3::new(
                    Fq::from(7592178415587923424u64),
                    Fq::from(9399622108558741788u64),
                    Fq::from(12654893332396084035u64),
                ),
                Fq3::new(
                    Fq::from(12678437828037174810u64),
                    Fq::from(2167899491043242876u64),
                    Fq::from(6190865625673193579u64),
                ),
                Fq3::new(
                    Fq::from(1480223285435494673u64),
                    Fq::from(17445458918271826306u64),
                    Fq::from(3497654872269564713u64),
                ),
                Fq3::new(
                    Fq::from(4394971171601069860u64),
                    Fq::from(16163573519134799361u64),
                    Fq::from(7407970631416741128u64),
                ),
                Fq3::new(
                    Fq::from(3904693916721222413u64),
                    Fq::from(14111409052569001470u64),
                    Fq::from(3944823175966235420u64),
                ),
                Fq3::new(
                    Fq::from(13361057081790561757u64),
                    Fq::from(14971804364859912820u64),
                    Fq::from(8381570998953724854u64),
                ),
                Fq3::new(
                    Fq::from(9610853259625319719u64),
                    Fq::from(9048842670925237621u64),
                    Fq::from(17507264183140614452u64),
                ),
                Fq3::new(
                    Fq::from(4711611854751477908u64),
                    Fq::from(5996535984055746175u64),
                    Fq::from(6093434598804609891u64),
                ),
            ]),
        ];

        assert_eq!(res, expected);
    }
}

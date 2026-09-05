# Notation rules for `algebra.jl/notebooks`

Status: local normative rule  
Scope: `lib/algebra.jl/notebooks/` の数学 notebook  
Owner context: 東京科学大学 鈴木咲衣研究室

この文書は、この directory 内の notebook にだけ適用する。workspace 全体の
記号規則ではない。

## 1. 流儀は冒頭で固定し、数式へ繰り返さない

R派・L派、left/right action、coproduct order などの流儀は notebook 冒頭で
文章と定義式により固定する。固定後の本文では、既定流儀を示すためだけの
添字を付けない。

例えば R派を既定値とした notebook では

$$
\Delta(E)=1\otimes E+E\otimes K,
\qquad
\Delta(F)=F\otimes1+K^{-1}\otimes F
$$

と書き、余積へ既定流儀の添字は付けない。同様に

$$
(E^iF^jH^k)_q^\vee,
\qquad E_q^\vee,
\qquad F_q^\vee,
\qquad H_q^\vee
$$

と書き、$(E^iF^jH^k)_{q,R}^\vee$ や $E_{q,R}^\vee$ とは書かない。

## 2. 複数流儀を同時比較するときだけラベルを付ける

既定流儀と別流儀を同じ箇所で比較する場合、既定側は無印のままとし、別側に
上付きの upright label を付ける。

$$
\Delta(E)=1\otimes E+E\otimes K,
\qquad
\Delta^{\mathrm L}(E)=E\otimes1+K\otimes E.
$$

流儀名を添字として各式へ常時伝播させない。

## 3. dual の記号

PBW basis に対する dual は一貫して $({\cdot})^\vee$ と書く。量子変形である
ことを区別する必要がある場合だけ $q$ を添字に置く。

- 古典: $(E^iF^jH^k)^\vee$
- 量子: $(E^iF^jH^k)_q^\vee$

同じ節で古典版との混同がなければ、量子版でも $q$ を省略してよい。ただし
一つの定理・表・計算ブロック内では表記を統一する。

## 4. coproduct の作用先

包絡代数側と座標環側の coproduct を同時に使うときは、流儀の R/L ではなく
作用先で区別する。

$$
\Delta_U:U\to U\otimes U,
\qquad
\Delta_O:O\to O\otimes O.
$$

作用先が文脈から明白なら単に $\Delta$ と書く。

## 5. 実装識別子

数学上省略した流儀ラベルを、理由なく関数名・変数名へ復活させない。複数流儀
を同じ program が実装し、dispatch または cache isolation に必要な場合だけ
`convention` を明示的な値として持たせる。

## 6. 適用時の確認

notebook を追加・修正するときは次を確認する。

1. 流儀が最初の定義箇所で一意に固定されている。
2. 既定流儀を示すだけの R/L 添字が残っていない。
3. R/L比較箇所以外で流儀ラベルを反復していない。
4. $\Delta_U$ と $\Delta_O$ の区別が必要な箇所で曖昧になっていない。
5. 同じ定理・表・code output 内で dual notation が統一されている。

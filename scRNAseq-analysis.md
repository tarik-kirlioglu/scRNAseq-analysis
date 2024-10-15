# scRNA-seq Analiz Projesi

## İçindekiler
- [Giriş](#giriş)
- [Bölüm 1](#bölüm-1)
  - [1.1. Kütüphaneleri Ortama Entegre Etme](#11-kütüphaneleri-ortama-entegre-etme)
  - [1.2. Veri Yükleme ve Seurat Nesnesi Oluşturma](#12-veri-yükleme-ve-seurat-nesnesi-oluşturma)
  - [1.3. Kalite Kontrol](#13-kalite-kontrol)
  - [1.4. Normalizasyon](#14-normalizasyon)
  - [1.5. Değişken Özelliklerin Belirlenmesi (özellik seçimi)](#15-değişken-özelliklerin-belirlenmesi-özellik-seçimi)
  - [1.6. Veri Ölçeklendirme](#16-veri-ölçeklendirme)
  - [1.7. Temel Bileşen Analizi ile Boyut Azaltma](#17-temel-bileşen-analizi-ile-boyut-azaltma)
  - [1.8. Hücreleri Kümeleme](#18-hücreleri-kümeleme)
  - [1.9. Non-linear Dimensional Reduction (UMAP/tSNE)](#19-non-linear-dimensional-reduction-umaptsne)
  - [1.10. Diferansiyel Olarak İfade Edilen Özelliklerin (Belirteçlerin) Belirlenmesi ve Hücre Tiplerini Atama](#110-diferansiyel-olarak-ifade-edilen-özelliklerin-belirteçlerin-belirlenmesi-ve-hücre-tiplerini-atama)
  - [1.11. Sonuçları Kaydetme](#111-sonuçları-kaydetme)


## Giriş

Tek hücreli RNA sekanslama (scRNA-seq), bireysel hücrelerin transkriptomlarını analiz etmek için kullanılan güçlü bir tekniktir. Bu yöntem, geleneksel bulk RNA-seq'in aksine, her bir hücrenin gen ekspresyon profilini ayrı ayrı incelememize olanak tanır. scRNA-seq, hücresel heterojenliği ortaya çıkarmak, nadir hücre tiplerini tanımlamak, hücre gelişim yollarını anlamak ve hastalık durumlarında hücresel değişiklikleri incelemek için kullanılır.

Şimdi, kodları ve açıklamalarını adım adım inceleyelim:

## Bölüm 1

## 1.1. Kütüphaneleri Ortama Entegre Etme

scRNAseq analizinde `Seurat` paketi ortama entegre edilmelidir. Aynı zamanda veri manipülasyonu ve görselleştirme için `tidyverse` paketi de faydalı olacaktır.

```R
library(Seurat)
library(tidyverse)
```

Paketlerin mevcut olmadığı durumda öncelikle `Seurat` ve `tidyverse` paketleri yüklenmelidir.

```R
install.packages('Seurat')
install.packages('tidyverse')
```

## 1.2. Veri Yükleme ve Seurat Nesnesi Oluşturma

Seurat, analiz boyunca verilerin saklanması ve işlenmesinde Seurat nesnesi adındaki veri türünü kullanır. Bu nedenle, ilk adım veriyi okumak ve Seurat nesnesi oluşturmaktır. Veri okuma adımında çeşitli okuma yöntemleri bulunmaktadır. `Read10x` fonksiyonu, 10x Genomics platformu tarafından üretilen çıktıları okumak için özel olarak tasarlanmıştır. Aynı zamanda, `ReadMtx` fonksiyonu kullanılarak da veriler okunur. `ReadMtx`, Daha esnektir. 10x verilerini okuyabilir, ancak diğer kaynaklardan gelen benzer formattaki verileri de okuyabilir.

Bu analizde, scRNA-seq veri setini okumak için Seurat paketinin `ReadMtx` fonksiyonu kullanılmıştır. `ReadMtx`, üç ayrı dosyadan oluşan bir veri setini okur:
`mtx`: Gene ekspresyon matrisini içeren dosya. Bu, hangi genin hangi hücrede ne kadar eksprese edildiğini gösteren sparse matristir.
`features`: Gen (veya diğer özelliklerin) isimlerini içeren dosya. Bu, matristeki satırların hangi genlere karşılık geldiğini belirtir.
`cells`: Hücre barkodlarını içeren dosya. Bu, matristeki sütunların hangi hücrelere karşılık geldiğini belirtir.

```R
counts <- ReadMtx(mtx = "GSE264489_RAW/GSM8218884_Dn1_matrix.mtx.gz",
               features = "GSE264489_RAW/GSM8218884_Dn1_features.tsv.gz",
               cells = "GSE264489_RAW/GSM8218884_Dn1_barcodes.tsv.gz")
```
Seurat nesnesi, `CreateSeuratObject` fonksiyonu ile oluşturulur. Bu nesne, tüm analiz boyunca kullanılacak ana veri yapısıdır. `min.cells` ve `min.features` parametreleri ayarlandığında verilere ilk filtreleme uygulanır. 

```R
seurat_obj <- CreateSeuratObject(counts = counts, project = "DN", min.cells = 3, min.features = 200)
```

## 1.3. Kalite Kontrol

Seurat nesnesi oluşturulduktan sonra kalite kontrol aşaması gelmektedir. Seurat, tek hücreli RNA sekanslama (scRNA-seq) verilerinin kalite kontrolü için çeşitli metrikler sunar. Bu metrikler, veri setinizdeki hücrelerin kalitesini değerlendirmenize ve düşük kaliteli hücreleri filtrelemenize olanak tanır. Seurat'ta sıkça kullanılan bazı önemli kalite kontrol metrikleri:

1.Hücre Başına Benzersiz Gen Sayısı (nFeature_RNA)
Hücrelerde kaç farklı genin ifade edildiğini belirten bu metrik, her hücredeki biyolojik çeşitliliği anlamak için kullanılır. Çok düşük gen sayısı, boş damlacıklar ya da düşük kaliteli hücrelere işaret edebilir. Öte yandan, aşırı yüksek bir gen sayısı, iki veya daha fazla hücrenin bir araya geldiği doublet gibi durumları gösterebilir.

2.Hücre Başına Toplam Molekül Sayısı (nCount_RNA)
Bu metrik, her hücredeki toplam RNA molekülü miktarını ölçer ve genellikle benzersiz gen sayısıyla doğrudan bir ilişkiye sahiptir, ancak hücre tipine bağlı olarak değişebilir. Güçlü bir korelasyon olması genellikle sağlıklı bir veri setine işaret eder, ancak bu korelasyon düşükse veri kalitesi sorgulanabilir. Düşük toplam molekül sayısı, boş damlacıklar veya kalitesiz hücreler ile ilişkilendirilebilirken, aşırı yüksek sayılar hücre doublet veya multipletlerin varlığını gösterebilir.

3.Mitokondriyal Gen Ekspresyon Oranı (percent.mt)
Mitokondriyal genlerin toplam gen ekspresyonuna oranını gösteren bu metrik, hücre stresini veya ölmekte olan hücreleri tespit etmek için kullanılır. Yüksek mitokondriyal gen ekspresyonu genellikle bu tür sorunlara işaret eder. 

Mitokondriyal genler genellikle "MT-" ile başlayan genlerdir ve bu genlerin yüzdesi, `PercentageFeatureSet` fonksiyonu kullanılarak hesaplanabilir.

```R
seurat_obj$percent.mt <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
```

Mitokondriyal gen ekspresyon oranı hesaplandıktan sonra, kalite kontrol metriklerini görselleştirmek için violin grafikleri kullanılabilir. Filtreleme için katı bir kural bulunmamakla birlikte, farklı veri setlerine göre bu adım değişkenlik gösterebilir. Genellikle violin grafikleri üzerinden hücre metriklerinin dağılımına bakarak, aşırı yüksek veya düşük değerlere sahip hücreler filtrelenir.

```R
VlnPlot(seurat_obj, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
```
![vlnplot_qc1](images/vlnplot_qc1.png)

Özellikler arasındaki ilişkileri incelemek için FeatureScatter kullanılabilir.
```R
plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
```

Bu veri seti için 400 ile 6000 arasında tespit edilen gen sayısı ve %10'dan düşük bir mitokondriyal ekspresyon yüzdesi kullanılmıştır. Fakat, farklı filtreleme eşikleri de kullanılabilir.

```R
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 400 & nFeature_RNA < 6000 & percent.mt < 10)
```
Filtreleme adımından sonra violin grafiği ile görselleştirme yapılabilir.

```R
VlnPlot(seurat_obj, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
```
![vlnplot_qc2](images/vlnplot_qc2.png)

Bazı durumlarda, temel filtreleme metriklerinin yanı sıra ek filtreleme adımlarına ihtiyaç duyulabilir. Özellikle en büyük zorluklardan biri, doublet ve multiplet hücrelerin varlığıdır. Standart filtreleme işlemleri sonucunda bazı doubletler veya multipletler fark edilmeyebilir. Bu tür hücreler veri setinde kaldığında analiz sonuçlarını yanıltabilir. Bu nedenle, doublet tespitine yönelik özel araçlar kullanmak önemlidir. 

[DoubletFinder](https://www.sciencedirect.com/science/article/pii/S2405471219300730), scRNA-seq veri setlerinde doublet hücreleri tespit etmek için tasarlanmış etkili bir araçtır. Bu yazılım, standart filtreleme işlemlerinden sonra bile fark edilmeyen problemli hücreleri belirleyerek veri setinin kalitesini artırır. DoubletFinder, yapay doubletler oluşturma ve karşılaştırma yöntemini kullanarak yüksek hassasiyetle çalışır. Bu sayede, analizlerde yanlış sonuçlara yol açabilecek doublet ve multipletlerin tespit edilmesine ve çıkarılmasına yardımcı olur.

## 1.4. Normalizasyon 

scRNA-seq (tek hücreli RNA dizileme) normalizasyonu, her hücredeki gen ekspresyon seviyelerinin karşılaştırılabilirliğini sağlamak için kritik bir adımdır. Bu aşama, farklı hücreler ve deneyler arasındaki teknik varyasyonları düzeltmek amacıyla yapılır, aynı zamanda biyolojik farklılıkları korumayı hedefler.

Normalizasyon sürecinde, her hücredeki toplam gen ekspresyonu, belirli bir referans değere (örneğin, toplam RNA miktarına) oranlanır. Yaygın bir yöntem olan CPM (Counts Per Million), her hücrenin gen ekspresyon seviyelerini, o hücredeki toplam gen ekspresyonuna böler ve bu oranları belirli bir çarpanla (genellikle 10,000 veya 1,000,000) çarpar. Bu yaklaşım, gen ekspresyon seviyelerini standart hale getirerek, farklı hücreler arasında biyolojik farklılıkların daha iyi değerlendirilmesini sağlar.

Normalizasyon işlemi, hücre boyutu ve RNA içeriği gibi faktörleri de dikkate alır. Ayrıca, genellikle normalizasyondan sonra logaritmik dönüşüm uygulanır (log-normalizasyon), bu da verilerin dağılımını iyileştirir ve downstream analizler için daha uygun hale getirir. Sonuç olarak, normalizasyon teknik varyasyonları azaltırken biyolojik farklılıkları koruyarak, scRNA-seq verilerinin doğru yorumlanmasını sağlar.

Normalizasyon için kullanılan `NormalizeData` fonksiyonunun çeşitli parametreleri bulunmakla birlikte, varsayılan ayarları çoğu durumda yeterli ve tercih edilir.

```R
seurat_obj <- NormalizeData(seurat_obj)
```


## 1.5 Değişken Özelliklerin Belirlenmesi (özellik seçimi)

Çok değişken özelliklerin belirlenmesi, scRNA-seq analizinde önemli bir adımdır ve genellikle özellik seçimi olarak adlandırılır. Bu süreç, veri setindeki genler arasındaki varyasyonu değerlendirmeye odaklanır ve en yüksek varyansa sahip genlerin belirlenmesini sağlar. Bu genler, biyolojik olarak anlamlı değişimleri yansıtma potansiyeline sahiptir ve analizde dikkate alınması gereken önemli özelliklerdir.

`FindVariableFeatures` fonksiyonu, scRNA-seq verilerinde en yüksek varyansa sahip genleri belirlemek için kullanılır. `selection.method` parametresi, en yüksek varyansa sahip genleri belirlemek için kullanılan yöntemi tanımlar ve varsayılan olarak "vst" (variance stabilizing transformation) olarak ayarlanmıştır; `nfeatures` parametresi ise seçilecek değişken özellik sayısını belirler ve varsayılan değeri 2000'dir. Genel olarak, 2000 gibi bir başlangıç değeri yaygın olarak kullanılır, ancak bu değeri veri setinin büyüklüğü, analiz hedefleri, genlerin varyans dağılımı ve hesaplama gücü gibi faktörlere göre özelleştirmek mümkündür.

```R
seurat_obj <- FindVariableFeatures(seurat_obj)
```

Sonuçlar, `VariableFeaturePlot` fonksiyonu ile görselleştirilebilir.

```R
top10 <- head(VariableFeatures(seurat_obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

```
![variable_features](images/variable_features.png)

## 1.6. Veri Ölçeklendirme 

scRNA-seq analizinde scale data (veri ölçeklendirme) aşaması, normalizasyondan sonra gelen bir adımdır. Bu süreç, seçilen özelliklerin veya genlerin ekspresyon değerlerini karşılaştırılabilir hale getirmeyi ve aşırı değişken genlerin etkisini dengelemeyi amaçlar. Böylece, yüksek ifade edilen genlerin baskınlığı önlenir. Ölçeklendirme işlemi, tipik olarak her genin ortalama ekspresyonunu 0'a, varyansını ise 1'e ayarlayarak gerçekleştirilir. Bu, genellikle z-skor dönüşümü kullanılarak yapılır.

```R
seurat_obj <- ScaleData(seurat_obj)
```
`ScaleData` fonksiyonu, verileri ölçeklendirirken aynı zamanda istenmeyen varyasyon kaynaklarını regrese etme (yani etkisini azaltma) imkanı sunar. Buradaki `vars.to.regress = "percent.mt"` parametresi, mitokondriyal gen yüzdesi ("percent.mt") ile ilişkili varyasyonun regrese edilmesini sağlar. Regresyon işlemi, teknik varyasyonu azaltmak için etkili olsa da, verilerin doğal biyolojik sinyallerini incelemek için analiz sürecinin başında yapılmaması daha uygundur, veriyi anladıktan sonra ihtiyaç duyulursa uygulanmalıdır.

```R
seurat_obj <- ScaleData(seurat_obj, vars.to.regress = "percent.mt")
```

Bu analizde geleneksel normalizasyon, değişken genlerin belirlenmesi ve veri ölçekleme yöntemleri kullanılmıştır. Bu işlemler sırasıyla `NormalizeData`, `FindVariableFeatures` ve `ScaleData` fonksiyonları ile gerçekleştirilmiştir. Alternatif olarak, [SCTransform](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1) yöntemi kullanılabilir; bu yöntem, normalizasyon, değişken genlerin belirlenmesi ve veri ölçekleme işlemlerini tek bir adımda gerçekleştiren daha gelişmiş bir yaklaşımdır. SCTransform, büyük ve heterojen veri setlerinde etkili sonuçlar verir; ancak bu yöntem, daha fazla hesaplama gücü ve zaman gerektirmesiyle birlikte karmaşıklığı nedeniyle zorlayıcı olabilir.

## 1.7. Temel Bileşen Analizi ile Boyut Azaltma

Boyut indirgeme, scRNA-seq verilerinde çok sayıda genin ölçüldüğü yüksek boyutlu veri setlerinin daha basit bir yapıya indirgenmesini sağlar. Yüksek boyutlu veriler, analiz edilmesi zor ve gürültüye duyarlı olabilir; bu nedenle boyut indirgeme ile verinin en önemli varyasyonlarını koruyarak daha az bileşenle temsil edilmesi gereklidir. Bu sayede analiz süreçleri hızlanır, görselleştirme kolaylaşır ve daha net biyolojik anlamlar çıkarılabilir.

Bu amaçla scRNA-seq analizinde en sık kullanılan yöntem Temel Bileşen Analizidir (PCA). PCA, verideki gen ekspresyonu varyasyonunu en iyi şekilde açıklayan "temel bileşenler" (principal components) adı verilen yeni eksenler oluşturur. Bu bileşenler, verideki en büyük varyansı açıklayan yönleri belirler ve veri, bu bileşenler üzerinden yeniden düzenlenir. `RunPCA` fonksiyonunda varsayılan olarak 50 temel bileşen (PCs) hesaplanır, ancak `npcs` parametresi ile bu sayı değiştirilebilir.

```R
seurat_obj <- RunPCA(seurat_obj)
```
PCA sonucunda, ilk birkaç temel bileşen verideki en fazla varyasyonu içerir. Genellikle, ilk bileşenler hücre tipleri arasındaki farklılıkları açıkça gösterirken, sonraki bileşenler daha küçük veya teknik varyasyonları yakalayabilir. Bu adım sayesinde, verideki biyolojik anlam taşıyan varyasyonlar ön plana çıkarılarak hücre kümeleri arasında ayrışma sağlanır.

PCA'dan sonra, kaç temel bileşenin (principal component) analiz için kullanılacağına karar vermek amacıyla `ElbowPlot` fonksiyonu kullanılır. Elbow plot, her bir bileşenin verideki toplam varyansı ne kadar açıkladığını görselleştirir. Grafik, bileşen sayısına karşı açıklanan varyansı gösterir ve genellikle bir "dirsek" (elbow) noktası fark edilir. Bu nokta, eklenen her yeni bileşenin varyansa katkısının azalmaya başladığı yerdir. Genellikle bu dirsek noktasında durmak, analiz için yeterli biyolojik bilgiyi içerir, ancak gereksiz bileşenleri eklemeden veri boyutunu da azaltır. Bu sayede, analizde en anlamlı ve varyansı yüksek bileşenler kullanılır.

```R
ElbowPlot(seurat_obj, ndims = 50)
```

![elbow_plot](images/elbow_plot.png)

Ayrıca, `DimHeatmap`, farklı temel bileşenlerin verideki hangi genetik varyasyonları temsil ettiğini ve varyasyonun bileşenler arasında nasıl dağıldığını incelemek için kullanılır. Bu analiz, kaç temel bileşenin seçileceğine karar verme sürecinde faydalı olabilir.

```R
DimHeatmap(seurat_obj, dims = 1:20, cells = 500, balanced = TRUE)
```

![heatmap_PCs](images/heatmap_PCs.png)

Elbow plot ve DimHeatmap sonuçlarına dayanarak analizde ilk 20 temel bileşen (PCs) seçildi. Kaç PCs seçileceğine karar vermek için kesin bir kural olmamakla birlikte, farklı sayıda PCs ile denemeler yapılarak en uygun sonuçlar belirlenebilir.

## 1.8. Hücreleri Kümeleme 
Hücreleri kümeleme adımı, scRNA-seq analizinde hücrelerin gen ekspresyon profillerine göre benzer gruplar halinde sınıflandırılmasını sağlar. Bu süreç, farklı hücre tiplerini, hücresel durumları veya biyolojik işlevleri ortaya çıkarmak için kullanılır. Genellikle, PCA gibi boyut indirgeme yöntemleri ile hücreler arasındaki temel farklılıklar belirlenir ve daha sonra bu farklılıklar baz alınarak hücreler kümelere ayrılır. Kümelenen hücreler, daha sonra biyolojik anlam taşıyan gruplara atanarak hücre tiplerinin tanımlanması ve yeni hücre popülasyonlarının keşfi sağlanır.

Kümeleme sürecinde ilk adım, hücreler arasındaki benzerliklerin belirlenmesi için `FindNeighbors` fonksiyonunun kullanılmasıdır. Bu fonksiyon, PCA ile elde edilen temel bileşenleri kullanarak her hücre arasındaki mesafeyi hesaplar ve hücreler arası ilişkiler üzerinden bir komşuluk grafiği oluşturur. Bu adım, hücrelerin birbirlerine ne kadar benzediğini belirlemek ve kümeleme için bir temel oluşturmak açısından kritiktir.

```R
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
```
Ardından, hücrelerin kümelere ayrılması için `FindClusters` fonksiyonu kullanılır. Bu fonksiyon, hücreler arası mesafelere dayalı olarak Louvain veya Leiden algoritmaları gibi kümeleme algoritmalarıyla hücreleri gruplar. Küme sayısı, `resolution` parametresi ile ayarlanabilir ve 0.1 ile 2.0 arasında değerler alabilir. Daha yüksek değerler daha fazla sayıda küçük kümeler oluştururken, daha düşük değerler daha az sayıda, daha büyük kümeler oluşturur.

```R
seurat_obj <- FindClusters(seurat_obj)
```

## 1.9. Non-linear Dimensional Reduction (UMAP/tSNE)

Non-linear boyut azaltma yöntemleri, özellikle karmaşık ve yüksek boyutlu verileri daha anlaşılır hale getirmek için kullanılır. scRNA-seq analizinde, UMAP (Uniform Manifold Approximation and Projection) ve t-SNE (t-distributed Stochastic Neighbor Embedding) gibi yöntemler, hücreler arasındaki benzerlikleri ve farklılıkları görselleştirmeye yardımcı olur.

```R
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
seurat_obj <- RunTSNE(seurat_obj, dims = 1:20)
```

Sonuçları görselleştirmek ve hücre kümelerini incelemek için `DimPlot` kullanılabilir. Bu fonksiyon, UMAP ve TSNE ile elde edilen düşük boyutlu temsilde hücrelerin dağılımını gösterir, böylece hücre kümeleri arasındaki ilişkileri ve farklılıkları daha iyi anlamak mümkün olur.

```R
dimplot1 <- DimPlot(seurat_obj, reduction = "umap")
dimplot1 <- DimPlot(seurat_obj, reduction = "tsne")

dimplot1 + dimplot2
```
![umap_vs_tsne](images/umap_and_tsne.png)

## 1.10. Diferansiyel Olarak İfade Edilen Özelliklerin (Belirteçlerin) Belirlenmesi ve Hücre Tiplerini Atama

Hücreleri kümelemek, scRNA-seq analizi sürecinin önemli bir aşamasıdır. Bu adım, benzer gen ekspresyon profillerine sahip hücrelerin gruplandırılmasını sağlar ve hücresel heterojenitenin anlaşılmasına katkıda bulunur. Kümeleme, hücrelerin biyolojik özelliklerine göre organize edilmesine olanak tanır ve genellikle belirteçler, literatür verileri veya otomatik yaklaşımlar kullanılarak gerçekleştirilir. Belirteçlere dayalı yöntemler, daha önce tanımlanmış gen setleri üzerinden hücre tiplerini belirlemeye yönelik olarak kullanılabilirken; otomatik yöntemler, makine öğrenimi algoritmaları ile hücreleri daha objektif bir şekilde sınıflandırma imkanı sunar.

Kümeleme adımının ardından, her bir küme için belirteçler atanarak bu hücre gruplarının tanımlanması sağlanır. Belirteçler, belirli hücre tiplerine özgü genlerin ekspresyon seviyeleri üzerinden belirlenir ve genellikle literatürdeki referanslar doğrultusunda seçilir. Belirteçlerden bazıları şunları içerir:

- CD4+ T cells: IL7R, CCR7
- CD8+ T cells: CD8A, CD8B
- CD14+ Monocytes: CD14, LYZ
- B cells: MS4A1, CD79A
- Natural killer (NK) cells: GNLY, NKG7, KLRB1
- Dendritic Cells: FCER1A, CST3
- FCGR3A+ Monocytes: FCGR3A, MS4A7
- Platelets: PPBP, PF4

Bu belirteçler, genel hücre tiplerinin ataması için bir referans sağlarken, alt hücre grupları için atamaların yapılabilmesi amacıyla ek belirteçlere ihtiyaç duyulabilir. Örneğin, T hücreleri CD4+ ve CD8+ alt gruplarına ayrılabilirken, B hücreleri de plazma hücreleri ve hafıza B hücreleri olarak sınıflandırılabilir. Böylece, hücresel heterojenite daha iyi anlaşılabilir ve spesifik işlevlerin incelenmesine olanak tanır.

Bu aşamada `FeaturePlot`, `DimPlot` ve `DoHeatmap` gibi görselleştirme araçları kullanılarak, belirteçlerin hücre grupları içindeki dağılımı ve yoğunluğu incelenir.

```R
FeaturePlot(seurat_obj, features = c("IL7R", "CCR7","CD8A", "CD8B", "CD14", "LYZ", "MS4A1", "CD79A", "GNLY",
                                     "NKG7", "FCER1A", "CST3", "PPBP", "PF4", "FCGR3A", "MS4A7"))

```
![featureplot_markers](images/featureplot_markers.png)

```R
DoHeatmap(seurat_obj, features = c("IL7R", "CCR7","CD8A", "CD8B", "CD14", "LYZ", "MS4A1", "CD79A", "GNLY",
                                   "NKG7", "FCER1A", "CST3", "PPBP", "PF4", "FCGR3A", "MS4A7")) + NoLegend()

```
![doheatmap_markers](images/doheatmap_markers.png)

Aynı zamanda, belirteçleri violin ve feature plot ile kombinleyerek ayrıntılı bir bakış sağlanabilir. Örneğin, CD8A ve CD8B belirteçleri küme 5'in güçlü belirteçleri olarak karşımıza çıkıyor.

![vlnplot_and_featureplot](images/vlnplot_and_featureplot.png)

Son olarak, diferansiyel olarak ifade edilen genlerin belirlenmesi için `FindAllMarkers` fonksiyonu kullanılır. Bu fonksiyon, her bir hücre kümesi için belirgin şekilde ifade edilen genleri tespit eder ve farklı hücre tipleri arasındaki gen ekspresyon farklılıklarını ortaya koyar. scRNA-seq verilerinin boyutunu göz önünde bulundurarak, kümelerdeki genlerin tespit oranını belirlemek için `min.pct` ve kat değişimini `logfc.threshold` yansıtmak için çeşitli parametreler girilebilir. Burada, her genin en az %25'inin hücrelerde ifade edilmesini sağlamak amacıyla `min.pct = 0.25` eşiğini belirler; ayrıca, genlerin en az %20'lik bir artış göstermesi gerektiğini belirtmek için `logfc.threshold = log(1.2)` parametresini kullanır.

```R
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))
```

Her küme için ifade edilen en yüksek on belirteç belirlenip, `DoHeatmap` fonksiyonu kullanılarak ısı haritaları oluşturulur. Bu görselleştirme, gen ekspresyonundaki farklılıkları ve kümeler arasındaki ilişkileri açık bir şekilde sunarak sonuçların yorumlanmasına yardımcı olur. Böylece, hücre tipleri ve bu tiplerin belirteçleri hakkında derinlemesine bir anlayış sağlanmış olur.

```R
top10_markers <- markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

DoHeatmap(seurat_obj, features = top10_markers$gene) + NoLegend() + theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 6))

```

![doheatmap_top10_markers](images/doheatmap_top10_markers.png)

CD4+ T hücrelerini temsil eden kümeler, belirteçler ve gen ekspresyon profillerine dayanarak küme 0, 3 ve 7 olarak belirlenmiştir. CD8A ve CD8B genlerinin güçlü ekspresyon gösterdiği küme 5 ise CD8+ T hücrelerine aittir. Küme 1, 2 ve 11 CD14+ monositlere, küme 9 ise FCGR3A+ monositlere atanmaktadır. MS4A1 ve CD79A genlerinin yüksek seviyede ifade edildiği küme 8 ve 13, B hücreleri olarak sınıflandırılmıştır. Küme 4 ve 6 NK hücrelerine, küme 10 ve 14 ise dendritik hücrelere karşılık gelmektedir. Son olarak, küme 15, plateletler olarak tespit edilmiştir.
 
Hücre kümelerinin bu şekilde atamaları yapılabilir; ancak bunun öznel olduğunu düşünebilirsiniz. Bu atamaları daha nesnel hale getirmek için SingleR gibi araçlar kullanılabilir. SingleR, referans veri setleri ile karşılaştırma yaparak hücre tiplerini otomatik olarak atamaya olanak tanır. Örneğin, referans olarak immün hücre veri tabanları kullanılarak, her bir hücrenin en yakın referans hücresi belirlenebilir ve bu sayede hücre tipi atamaları daha sistematik ve objektif bir şekilde gerçekleştirilebilir. SingleR çıktıları, manuel analizlerle karşılaştırıldığında tutarlılık sağlanması da analiz doğruluğunu artırır.

Hücre tipleri, hücre kümelerine atanabilir. Her küme için belirlenen hücre tipleri bir vektöre aktarılır ve names fonksiyonu ile bu hücre tipleri, Seurat objesindeki mevcut kümelerle eşleştirilir. `RenameIdents` fonksiyonu kullanılarak, Seurat objesindeki hücre kümeleri bu yeni hücre tiplerine göre yeniden adlandırılır. Bu sayede, her hücre kümesi biyolojik olarak anlamlı hücre tipleriyle eşleştirilmiş olur.

```R
new.ids <- c("CD4+ T", "CD14+ Mono", "CD14+ Mono", "CD4+ T", "NK", "CD8+ T", "NK", "CD4+ T", "B cell",
             "FCGR3A+ Mono", "DC", "CD14+ Mono", "NK", "B cell", "DC", "Platelets" )
names(new.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new.ids)
```

`DimPlot` fonksiyonu kullanılarak atanan hücre tipleriyle bir görselleştirme yapılabilir.

```R
DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 1) + NoLegend()
```
![dimplot_cell_type](images/dimplot_cell_type.png)

## 1.11. Sonuçları Kaydetme

Buraya kadar yapılan analizler, tek bir scRNA-seq verisi ile gerçekleştirilen analizleri kapsar. Bu işlemin sonunda, Seurat objesini kaydederek yeniden analiz yapmaya gerek kalmaz. Seurat objesi, `saveRDS` fonksiyonu kullanılarak kaydedilebilir ve daha sonra `readRDS` fonksiyonu ile tekrar yüklenebilir. Bu sayede, analizin ilerleyen aşamalarında aynı obje üzerinde çalışmaya devam edilebilir.

```R
saveRDS(seurat, file="DN1/seurat_obj_final.rds")
seurat <- readRDS("DN1/seurat_obj_final.rds")
```

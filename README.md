# article_13778_Polarimetric_Decomposition
日本発の衛星データプラットフォーム Tellus のオウンドメディア「宙畑」の記事、https://sorabatake.jp/13778 で利用しているコードです。
偏波の成分分解を用いて、SARの偏波画像から構造物と自然物を選り分けることが出来るのかを検証しました。

なお、SARデータの取得にはTellusの開発環境を申し込む必要があります。

## ソースコード(./src/配下を参照)
- 00.py
  - CEOSのダウンロード
- 01.py
  - CEOS→SLCに変換
- 02.py
  - パウリ分解
- 03.py
  - パウリ分解+freeman-duran
- slcinfo.py
  - SLCの構造体

## ライセンス、利用規約
ソースコードのライセンスは CC0-1.0（Creative Commons Zero v1.0 Universal）ライセンスです。  
今回コード内で PALSAR-2 データを用いております。利用ポリシーは以下をご参考下さい。
https://www.tellusxdp.com/market/tool_detail/de3c41ac-a8ca-4170-9028-c9e1a39841e1/e364c31c-bfad-49d0-bd6d-f2bc11d67386
※サイトの閲覧にはTellusへのログインが必要です。

## 貢献方法
プルリクエストや Issue はいつでも歓迎します。



by charmegiddo

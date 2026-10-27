# QTSwE 後端 (C)

依 Kuo, Chou & Chen, *Quantum-inspired algorithm for cyber-physical visual surveillance deployment systems*, Computer Networks 117 (2017) 的 QTSwE 邏輯，對任意 m×n 範圍自動求「滿足 Pg ≥ Tg（且可選：全部連回 Sink）的最少感測器擺法」。

## Windows（直接在 VS Code 的 PowerShell）
1. 安裝 MSYS2（https://www.msys2.org），開 "MSYS2 UCRT64" 執行 `pacman -S mingw-w64-ucrt-x86_64-gcc`，再把 `C:\msys64\ucrt64\bin` 加到 Windows 的 PATH。
2. 在 VS Code 終端機（PowerShell）進入此資料夾，執行 `.\build.bat`，成功會產生 `qtswe.exe`。
3. 啟動：`.\qtswe.exe serve`；命令列測試：`.\qtswe.exe run --rows 30 --cols 30 --tg 0.3`

## 編譯與執行（Mac / Linux）
```bash
make                      # 需要 gcc；預設用 OpenMP 多執行緒（make OPENMP=0 可關閉）
./qtswe run --rows 50 --cols 50 --tg 0.3              # 命令列，印出 ASCII 地圖
./qtswe run --rows 120 --cols 90 --tg 0.6 --json      # 輸出 JSON
./qtswe serve --port 8080                              # 啟動給網頁用的 HTTP 後端
./qtswe bench --rows 30 --cols 30 --tg 0.3 --runs 10   # 消融實驗，CSV 輸出到 stdout
python3 verify.py < result.json                         # 用 Python 獨立重算 Pg 與連通性
```

## HTTP API（CORS 已開）
`GET /health`　·　`GET /optimize?rows=25&cols=40&tg=0.3&conn=1&seed=1`
可用參數：`rows cols tg rcon rmax alpha beta rt conn sink_r sink_c(1起算) seed pop gens detect delta localsize lradius entangle local retries k(固定感測器數) time(秒) include_pg`
回傳 JSON：`sensors`（[列,行]，1 起算）、`k`、`feasible`、`covered`、`connected`、`pg_min`、`attempts`（每個試過的 k）、`curve`（適應度曲線）…

## 演算法對照論文
| 論文 | 本實作 |
|---|---|
| Q 矩陣 \|β\|²、Measure（R>\|β\|² 則為 1） | 同；每代 p=100 組 |
| 糾纏啟發初始化（沿軸依適當距離排鏈 + 隨機擾動） | 三角晶格鏈（交錯排列）+ 隨機擾動，間距由 k 決定；Tg 對應的理論間距 `s_prop` 用來估起始 k |
| 適應度 f=w1·f1+w2·f2（達標格數 + 連回 Sink 的感測器數），軟限制 | 同，固定 k；f1=RC 且 f2=k 才算可行 |
| 以最佳解吸引、最差解排斥更新 Q（旋轉 0.1） | 全域最佳 vs 當代最差，delta=0.1 |
| 停滯 10 代 → 區域搜尋 | 同；每顆感測器往鄰近格試移動 |
| 最少感測器數 | 外層：由晶格估計值往上找到可行，再每次減 1（以前一個可行解為暖啟動）直到修不回來 |

### 與論文不同／未做的地方（請在報告中如實說明）
- **沒有實作禁忌表 (tabu list)**，只有「停滯觸發的區域搜尋」。
- 固定感測器數 k 時，Q 初始值設為 1−k/(RC)（不是 0.5），Measure 後再修補成剛好 k 顆。
- 區域搜尋預設 `lradius=2`（5×5 範圍）並允許同分側移；論文是 8 鄰格，用 `--lradius 1` 可還原。
- 外層「逐顆減少」與暖啟動是本實作加的，論文沒有此流程。
- Sink 預設在地圖中心且不可放感測器（論文未指定）。
- 論文 Section 6 文字寫 0.3、Table 13 寫旋轉 0.1π，本實作採 0.1。
- 50×50, Rcon=1.414, Rmax=5, Rt=5 的實測（seed=1）：Tg=0.3/0.6/0.9 無連通 58/100/212 顆（論文 62/115/220）；Tg=0.3 含連通 72 顆（論文 83）。數字比論文少，可能是幾何細節（距離、邊界）不同，**不能宣稱完全重現論文**。

## 退火 (SA) 接受機制（可選，預設關閉）
- `--sa_ref 1`（A）：每代拿「當代最好」與「參考解」比，變差時以 exp(-delta/T) 機率接受，被接受的解會成為吸引 Q 的參考解。另外保存「歷史最好解」，最終答案與可行判斷都看歷史最好解。
- `--sa_ls 1`（B）：區域搜尋時，若某顆感測器所有鄰近挪法都比現在差，以 exp(-delta/T) 機率接受其中一個隨機挪法；整段區域搜尋回傳途中看過最好的擺法。
- 溫度每次嘗試從 `--sa_t0`（預設 2）幾何降到 `--sa_tend`（預設 0.1）。A、B 共用同一個溫度。
- 兩者都關閉時，結果與未加入此功能前逐位元相同。
- `bench` 會多跑 full+SA_A / full+SA_B / full+SA_AB 三組。

## 固定感測器數的成功率比較（消融用）
```
./qtswe bench --rows 50 --cols 50 --tg 0.3 --conn 0 --k 58 --runs 12
```
加上 `--k` 後，每組只測「只用 k 顆能不能全部達標」（從零開始，最多 gens 代），輸出每組的成功次數、平均代數、平均覆蓋格數；CSV 在 stdout，摘要表在 stderr。接近極限的 k 最容易看出差別。

## 給下一步（ADQISA）預留的位置
- `qp_t.on_attempt(k, feasible, gens_used, ud)`：每試完一個 k 就呼叫，回傳非 0 可提早停止整個搜尋 → ADQISA 的「還值不值得繼續降 k」判斷可放這裡。
- `qp_t.on_generation(k, gen, best_f, ud)`：每代回報 → 可做收斂／停滯判斷。
- 消融開關：`entangle`、`local`、`lradius`、`retries`、`delta`、`pop`、`gens`；`bench` 子指令已內建 full / no_entangle / no_local / plain_QTS 四組（同 seed 序列）。
- 大範圍：上限 4,000,000 格；固定 `--k` 可只測某個數量的成功率；`--time` 設總時間上限（HTTP 預設 120 秒）。


## Result cache

`qtswe serve` stores every finished `/optimize` result in `cache/<hash>.json` (next to where you start it).
Asking again with the same parameters (including `seed`) returns the stored layout in milliseconds
(`"cached":true`; `seconds` is the original compute time). Not cached: runs that hit the time limit, `seed=0`.
`threads` and `time` are not part of the key. Add `&nocache=1` to force a recompute; delete the `cache/` folder to clear it.
`qtswe serve --cache DIR` changes the folder, `--cache off` disables it. Delete the cache after rebuilding with
algorithm changes, otherwise old results are replayed.

# Zenodo PI 数据断点续传（修复作用域，路径作参数）
function Get-File($outdir, $name, $url) {
  $dest = Join-Path $outdir $name
  for ($i = 1; $i -le 3000; $i++) {
    & curl.exe -L -C - -s --connect-timeout 20 --max-time 900 -o $dest $url
    $c = $LASTEXITCODE
    if ($c -eq 0 -or $c -eq 33) { Write-Output "$name done at try $i"; break }
    Start-Sleep -Milliseconds 500
  }
}

$out = $PSScriptRoot
Get-File $out "PI_sc_h5.zip" "https://zenodo.org/api/records/19697597/files/filtered%20h5%20files%20for%20peri-implantitis%20scRNA-seq.zip/content"
Get-File $out "PI_spaceranger.zip" "https://zenodo.org/api/records/19697597/files/space%20ranger%20output.zip/content"
Write-Output "ZENODO_ALL_DONE"

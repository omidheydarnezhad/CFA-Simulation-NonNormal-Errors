# --- 0. نصب و بارگذاری پکیج‌ها ---
# این بخش را فقط یک بار در ابتدای پروژه اجرا کنید.
# install.packages(c("psych", "EGAnet", "tidyverse", "parallel", "MASS", "Matrix", "nFactors"), dependencies = TRUE)

# این بخش را در ابتدای هر سشن R اجرا کنید.
library(psych)
library(EGAnet)
library(tidyverse)
library(parallel)
library(MASS)    # برای mvrnorm
library(Matrix)  # برای nearPD
library(nFactors) # برای nScree

# --- 1. تعریف تابع generate_data_improved ---
# این تابع داده‌های مشاهده‌شده را با ساختار عاملی مشخص و نویز کنترل‌شده تولید می‌کند.
generate_data_improved <- function(n_items, n_factors, n_obs, loading_strength_range, factor_corr_val, noise_variance = 0.1, cross_loading_prob = 0.0) {
  
  # 1. تعریف ماتریس بارهای عاملی (Lambda) با ساختار ساده
  # هر آیتم فقط روی یک عامل بار معنادار دارد و بارهای متقاطع به صورت تصادفی اضافه می‌شوند.
  lambda <- matrix(0, nrow = n_items, ncol = n_factors)
  
  # توزیع آیتم‌ها بین عوامل
  items_per_factor_base <- floor(n_items / n_factors)
  remainder_items <- n_items %% n_factors
  
  item_idx <- 1 # اندیس شروع آیتم برای عامل فعلی
  for (f in 1:n_factors) {
    num_items_this_factor <- items_per_factor_base + (f <= remainder_items)
    
    # تعیین آیتم‌های مربوط به این عامل
    current_items_indices <- item_idx:(item_idx + num_items_this_factor - 1)
    
    # تولید بارهای عاملی اصلی
    lambda[current_items_indices, f] <- 
      runif(num_items_this_factor, min = loading_strength_range[1], max = loading_strength_range[2])
    
    # اضافه کردن بارهای متقاطع برای آیتم‌های این عامل
    for (j in current_items_indices) {
      for (other_f in setdiff(1:n_factors, f)) { # برای عوامل دیگر
        if (runif(1) < cross_loading_prob) {
          # بار متقاطع ضعیف (مثلا 0.1 تا 0.3)
          # اطمینان از اینکه بار متقاطع خیلی بالا نباشد
          lambda[j, other_f] <- runif(1, min = 0.1, max = min(0.3, loading_strength_range[1] - 0.05)) # اطمینان از پایین‌تر بودن از بار اصلی
        }
      }
    }
    item_idx <- item_idx + num_items_this_factor
  }
  
  # 2. تعریف ماتریس همبستگی عوامل (Phi)
  if (factor_corr_val == 0) {
    phi <- diag(n_factors) # Orthogonal (عوامل مستقل)
  } else {
    phi <- matrix(factor_corr_val, nrow = n_factors, ncol = n_factors)
    diag(phi) <- 1 # Oblique (عوامل همبسته)
  }
  
  # اطمینان از مثبت-معین بودن ماتریس همبستگی (ضروری برای mvrnorm)
  # این مرحله از خطاهای احتمالی در mvrnorm جلوگیری می‌کند.
  phi <- as.matrix(nearPD(phi)$mat) 
  
  # 3. تولید نمرات عامل (F)
  factors_scores <- mvrnorm(n = n_obs, mu = rep(0, n_factors), Sigma = phi)
  
  # 4. تولید واریانس خطای اندازه‌گیری (Uniqueness)
  # communality (اشتراک) برای هر آیتم: h2 = sum(lambda_ij^2) over j
  # u2_base = 1 - h2. این واریانس منحصر به فرد (خطای اندازه‌گیری بدون نویز اضافی) است.
  h2 <- rowSums(lambda^2)
  u2_base <- 1 - h2 
  
  # اضافه کردن نویز به واریانس منحصربه‌فرد
  # total uniqueness = u2_base + noise_variance
  u2 <- u2_base + noise_variance 
  
  # اطمینان از اینکه واریانس خطاها مثبت است
  # اگر h2 به 1 نزدیک باشد، u2_base نزدیک به صفر می‌شود.
  # اگر u2 منفی شود (h2 > 1)، آن را به حداقل مقدار مثبت محدود می‌کنیم.
  u2[u2 <= 0] <- 0.01 # حداقل واریانس خطا را تضمین می‌کند
  
  # 5. تولید خطاهای اندازه‌گیری (Epsilon)
  # خطاها از توزیع نرمال با میانگین صفر و انحراف معیار sqrt(u2)
  epsilon <- matrix(rnorm(n = n_obs * n_items, mean = 0, sd = rep(sqrt(u2), each = n_obs)),
                    nrow = n_obs, ncol = n_items, byrow = TRUE)
  
  # 6. تولید داده‌های مشاهده‌شده (X)
  data_sim <- factors_scores %*% t(lambda) + epsilon
  
  # نام‌گذاری ستون‌ها برای خوانایی بهتر
  colnames(data_sim) <- paste0("item_", 1:n_items)
  
  return(as.data.frame(data_sim))
}

# --- 2. تعریف تابع apply_methods ---
# این تابع داده‌های شبیه‌سازی‌شده را می‌گیرد و روش‌های مختلف حفظ عامل را روی آن اعمال می‌کند.
apply_methods <- function(data) {
  # نرمال‌سازی داده‌ها (ضروری برای برخی روش‌ها)
  data_scaled <- scale(data)
  
  # 1. Parallel Analysis (PA)
  # n.iter افزایش یافته برای پایداری بیشتر
  pa_result <- tryCatch({ 
    fa.parallel(data_scaled, fa = "fa", n.iter = 500, plot = FALSE)$nfact 
  }, error = function(e) NA)
  
  # 2. Velicer's MAP Test
  map_result <- tryCatch({ 
    VSS_output <- VSS(data_scaled, n = ncol(data_scaled), plot = FALSE, fm="pa")
    VSS_output$map[1] # Velicer's MAP test (اولین مقدار استاندارد)
  }, error = function(e) NA)
  
  # 3. Kaiser Criterion (Eigenvalue > 1)
  kaiser_result <- tryCatch({ 
    sum(eigen(cor(data_scaled))$values > 1) 
  }, error = function(e) NA)
  
  # 4. Exploratory Graph Analysis (EGA)
  ega_result <- tryCatch({ 
    EGA(data_scaled, model = "glasso", plot.EGA = FALSE)$n.dim 
  }, error = function(e) NA)

  # 5. Scree Plot (automatic interpretation via nFactors)
  scree_result <- tryCatch({
    ev <- eigen(cor(data_scaled))$values
    nScree(x = ev)$Components$nScree
  }, error = function(e) NA)
  
  # 6. BIC / AIC
  min_bic <- Inf
  min_aic <- Inf
  bic_suggested <- NA
  aic_suggested <- NA
  
  # حداکثر تعداد عوامل ممکن برای بررسی
  # حداقل 1 و حداکثر min(تعداد آیتم‌ها - 1، یک سوم حجم نمونه)
  max_factors_to_test <- min(ncol(data_scaled) - 1, floor(nrow(data_scaled) / 3)) 
  if (max_factors_to_test < 1) max_factors_to_test <- 1 # اطمینان از بررسی حداقل 1 عامل
  
  for (nf in 1:max_factors_to_test) {
    tryCatch({
      # اجرای تحلیل عاملی برای تعداد عوامل مشخص
      fa_model <- fa(data_scaled, nfactors = nf, rotate = "none", fm = "pa")
      
      # به‌روزرسانی BIC و AIC اگر مدل بهبود یابد
      if (!is.null(fa_model$BIC) && fa_model$BIC < min_bic) {
        min_bic <- fa_model$BIC
        bic_suggested <- nf
      }
      if (!is.null(fa_model$AIC) && fa_model$AIC < min_aic) {
        min_aic <- fa_model$AIC
        aic_suggested <- nf
      }
    }, error = function(e) {
      # در صورت بروز خطا در fa() برای تعداد عوامل خاص، آن را نادیده بگیرید
      # (این معمولاً برای جلوگیری از توقف شبیه‌سازی در موارد خاص است)
    })
  }

  return(list(
    pa = pa_result,
    map = map_result,
    kaiser = kaiser_result,
    ega = ega_result,
    scree = scree_result,
    bic = bic_suggested,
    aic = aic_suggested
  ))
}

# --- 3. تعریف تابع run_replication ---
# این تابع یک تکرار واحد از شبیه‌سازی را برای یک سناریوی خاص اجرا می‌کند.
run_replication <- function(rep_id, scenario) {
  # تولید داده
  data <- tryCatch({
    generate_data_improved(
      n_items = scenario$n_items,
      n_factors = scenario$n_factors_true, 
      n_obs = scenario$sample_size,
      loading_strength_range = unlist(scenario$loading_strength_range), # unlist برای دسترسی به مقادیر min/max
      factor_corr_val = scenario$factor_corr,
      noise_variance = scenario$noise_variance,
      cross_loading_prob = scenario$cross_loading_prob
    )
  }, error = function(e) { 
    warning(paste("خطا در تولید داده برای سناریو", paste(scenario, collapse="-"), "تکرار", rep_id, ":", e$message))
    NULL # در صورت خطا، NULL برگردان
  })
  
  if (is.null(data)) return(NULL) # اگر تولید داده ناموفق بود، NULL برگردانده و ادامه ندهد
  
  # اجرای روش‌های تعیین عامل
  methods_results <- apply_methods(data)
  
  # ساخت یک تیبل برای نتایج این تکرار
  result_row <- tibble(
    n_factors_true = scenario$n_factors_true,
    n_items = scenario$n_items,
    sample_size = scenario$sample_size,
    loading_min = unlist(scenario$loading_strength_range)[1], # ذخیره بازه بار عاملی
    loading_max = unlist(scenario$loading_strength_range)[2],
    factor_corr = scenario$factor_corr,
    noise_variance = scenario$noise_variance,
    cross_loading_prob = scenario$cross_loading_prob,
    replication = rep_id,
    pa = methods_results$pa,
    map = methods_results$map,
    kaiser = methods_results$kaiser,
    ega = methods_results$ega,
    scree = methods_results$scree,
    bic = methods_results$bic,
    aic = methods_results$aic
  )
  
  return(result_row)
}

# --- 4. تنظیم و اجرای شبیه‌سازی ---
# تنظیم seed برای تکرارپذیری
set.seed(123)

# تعریف Design Grid (سناریوها)
design_grid <- expand.grid(
  n_factors_true = c(2, 3, 5),
  n_items = c(15, 30, 50),
  sample_size = c(50, 100, 300, 1000), # شامل نمونه‌های کوچک‌تر
  # تعریف بازه‌های بارهای عاملی به عنوان لیست برای استفاده صحیح در expand.grid
  loading_strength_range = list(c(0.3, 0.5), c(0.5, 0.7), c(0.7, 0.9)), 
  factor_corr = c(0.0, 0.3, 0.5), # 0.0 برای متعامد، 0.3 و 0.5 برای مایل
  noise_variance = c(0.0, 0.1, 0.3), # 0.0 برای بدون نویز، 0.1 (کم)، 0.3 (متوسط)
  cross_loading_prob = c(0.0, 0.2, 0.4) # 0.0 برای بدون بارهای متقاطع، 0.2 (کم)، 0.4 (متوسط)
)

n_replications <- 1000 # تعداد تکرارها برای هر سناریو. (برای مقاله: حداقل 500-1000)

# تنظیمات موازی‌سازی
n_cores <- detectCores() - 1 # استفاده از تمام هسته‌ها به جز یکی
if (n_cores < 1) n_cores <- 1 # حداقل یک هسته را تضمین می‌کند
cl <- makeCluster(n_cores)

# اکسپورت کردن توابع و پکیج‌های لازم به کلاسترها
clusterExport(cl, c("generate_data_improved", "apply_methods", "run_replication", "n_replications"))
clusterEvalQ(cl, {
  # اطمینان از بارگذاری پکیج‌ها در هر هسته
  library(psych)
  library(EGAnet)
  library(tidyverse)
  library(MASS)
  library(Matrix)
  library(nFactors)
})

# اجرای شبیه‌سازی به‌صورت موازی
all_results_list <- parLapply(cl, 1:nrow(design_grid), function(i) {
  scenario <- design_grid[i, ]
  
  # چاپ پیشرفت در کنسول (با اطلاعات سناریو)
  cat(paste0("Running Scenario ", i, " of ", nrow(design_grid), 
             " (True Factors: ", scenario$n_factors_true, 
             ", Items: ", scenario$n_items, 
             ", N: ", scenario$sample_size, 
             ", Loadings: ", unlist(scenario$loading_strength_range)[1], "-", unlist(scenario$loading_strength_range)[2],
             ", Corr: ", scenario$factor_corr, 
             ", Noise: ", scenario$noise_variance,
             ", Cross-Load: ", scenario$cross_loading_prob, ")\n"))
  
  # اجرای تکرارها برای این سناریو
  scenario_results <- lapply(1:n_replications, function(rep) {
    run_replication(rep, scenario)
  })
  
  # ترکیب نتایج تکرارها برای این سناریو
  bind_rows(scenario_results)
})

# بستن خوشه پس از اتمام شبیه‌سازی
stopCluster(cl)

# ترکیب و ذخیره‌سازی تمامی نتایج نهایی
final_results_df <- bind_rows(all_results_list)

# ذخیره نتایج خام (برای تحلیل‌های بعدی و تکرارپذیری)
write.csv(final_results_df, "simulation_results_raw.csv", row.names = FALSE)
saveRDS(final_results_df, "simulation_results_raw.rds")

cat("شبیه‌سازی کامل شد و نتایج در simulation_results_raw.csv ذخیره شد.\n")

# --- 5. تحلیل و مصورسازی نتایج ---
# بارگذاری نتایج (اگر سشن جدیدی شروع کرده‌اید یا فایل را جابجا کرده‌اید)
# final_results_df <- readRDS("simulation_results_raw.rds")

# تعریف نام روش‌ها برای راحتی کار در تحلیل
methods_to_evaluate <- c("pa", "map", "kaiser", "ega", "scree", "bic", "aic")

# محاسبه دقت (Accuracy)، خطای مطلق میانگین (MAE) و سوگیری (Bias)
summary_metrics <- final_results_df %>%
  pivot_longer(
    cols = all_of(methods_to_evaluate),
    names_to = "method",
    values_to = "estimated_factors"
  ) %>%
  mutate(
    # محاسبه درست بودن تخمین (1 اگر صحیح باشد، 0 در غیر این صورت)
    correct = as.integer(estimated_factors == n_factors_true),
    # محاسبه خطای مطلق (|پیش‌بینی شده - واقعی|)
    abs_error = abs(estimated_factors - n_factors_true),
    # محاسبه سوگیری (پیش‌بینی شده - واقعی)
    bias = estimated_factors - n_factors_true
  ) %>%
  group_by(method, n_factors_true, n_items, sample_size, loading_min, loading_max, factor_corr, noise_variance, cross_loading_prob) %>%
  summarise(
    accuracy = mean(correct, na.rm = TRUE), # میانگین دقت
    mae = mean(abs_error, na.rm = TRUE),     # میانگین خطای مطلق
    mean_bias = mean(bias, na.rm = TRUE),    # میانگین سوگیری
    n_na_results = sum(is.na(estimated_factors)), # تعداد NA برای هر روش در هر سناریو
    .groups = "drop" # حذف گروه‌بندی
  )

# نمایش و ذخیره نتایج خلاصه شده
print(summary_metrics)
write.csv(summary_metrics, "simulation_summary_metrics.csv", row.names = FALSE)
saveRDS(summary_metrics, "simulation_summary_metrics.rds")

cat("تحلیل نتایج کامل شد و metrics در simulation_summary_metrics.csv ذخیره شد.\n")

# --- 6. رسم نمودارها (نمونه‌هایی برای شروع تحلیل) ---
# توجه: با تعداد زیاد سناریوها، بهتر است نمودارها را برای چند سناریوی کلیدی فیلتر کنید
# یا از فاکتوربندی‌های دقیق برای نمایش همه استفاده کنید.

# نمودار ۱: دقت بر حسب حجم نمونه، با تفکیک بر اساس تعداد عوامل واقعی و همبستگی عوامل
ggplot(summary_metrics, aes(x = sample_size, y = accuracy, color = method, group = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  # فیلتر کردن برای یک زیرمجموعه از سناریوها برای خوانایی بهتر نمودار
  filter(noise_variance == 0.1, cross_loading_prob == 0.0, n_items == 30, loading_min == 0.5) +
  facet_grid(n_factors_true ~ factor_corr, 
             labeller = labeller(n_factors_true = "True Factors: {.}", factor_corr = "Factor Corr: {.}")) +
  labs(
    title = "Accuracy of Factor Retention Methods by Sample Size (Selected Scenarios)",
    x = "Sample Size",
    y = "Accuracy Rate (Proportion Correct)",
    color = "Method"
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.text = element_text(size = 8, face = "bold")) +
  scale_x_continuous(breaks = unique(summary_metrics$sample_size)) +
  scale_color_brewer(palette = "Set1")

# نمودار ۲: MAE (Mean Absolute Error) بر حسب Noise Variance
ggplot(summary_metrics, aes(x = factor(noise_variance), y = mae, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  # فیلتر کردن برای یک زیرمجموعه
  filter(n_factors_true == 3, sample_size == 300, n_items == 30, factor_corr == 0.0, cross_loading_prob == 0.0) +
  labs(
    title = "Mean Absolute Error by Noise Variance (Selected Scenario)",
    x = "Noise Variance",
    y = "Mean Absolute Error",
    fill = "Method"
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set1")

# نمودار ۳: دقت بر حسب Cross-Loading Probability
ggplot(summary_metrics, aes(x = factor(cross_loading_prob), y = accuracy, color = method, group = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  # فیلتر کردن برای یک زیرمجموعه
  filter(n_factors_true == 3, sample_size == 300, n_items == 30, factor_corr == 0.0, noise_variance == 0.1) +
  labs(
    title = "Accuracy by Cross-Loading Probability (Selected Scenario)",
    x = "Cross-Loading Probability",
    y = "Accuracy Rate",
    color = "Method"
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_brewer(palette = "Set1")

# نمودار ۴: سوگیری (Bias) بر حسب حجم نمونه
ggplot(summary_metrics, aes(x = sample_size, y = mean_bias, color = method, group = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") + # خط مرجع برای سوگیری صفر
  filter(n_factors_true == 3, n_items == 30, loading_min == 0.5, factor_corr == 0.0, noise_variance == 0.1, cross_loading_prob == 0.0) +
  labs(
    title = "Bias of Factor Retention Methods by Sample Size (Selected Scenario)",
    x = "Sample Size",
    y = "Mean Bias (Estimated - True)",
    color = "Method"
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_x_continuous(breaks = unique(summary_metrics$sample_size)) +
  scale_color_brewer(palette = "Set1")

# نمودار ۵: دقت بر حسب بازه بارهای عاملی
ggplot(summary_metrics, aes(x = factor(loading_min), y = accuracy, color = method, group = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  filter(sample_size == 300, n_factors_true == 3, factor_corr == 0.0, noise_variance == 0.1, cross_loading_prob == 0.0) +
  labs(
    title = "Accuracy by Loading Strength Range (Selected Scenario)",
    x = "Minimum Loading Strength",
    y = "Accuracy Rate",
    color = "Method"
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_brewer(palette = "Set1")
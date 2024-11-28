// ==UserScript==
// @name         PubMed 页面显示期刊影响因子
// @namespace    http://tampermonkey.net/
// @version      0.6.2
// @description  在 PubMed 和 PMC 页面显示期刊影响因子
// @author       Your name
// @match        https://pubmed.ncbi.nlm.nih.gov/*
// @match        https://pmc.ncbi.nlm.nih.gov/articles/*
// @updateURL    https://raw.githubusercontent.com/ydx1013/DATA/refs/heads/master/pubmed.js
// @downloadURL  https://raw.githubusercontent.com/ydx1013/DATA/refs/heads/master/pubmed.js
// @grant        none
// ==/UserScript==

(function() {
    'use strict';
    console.log('油猴脚本已加载');

    // 预定义颜色池
    const colorSchemes = [
        { bg: '#e3f2fd', border: '#e3f2fd', text: '#1565C0' },  // 蓝色
        { bg: '#e0f7fa', border: '#e0f7fa', text: '#00838F' },  // 青色
        { bg: '#e8f5e9', border: '#e8f5e9', text: '#2E7D32' },  // 绿色
        { bg: '#f3e5f5', border: '#f3e5f5', text: '#6A1B9A' },  // 紫色
        { bg: '#fff3e0', border: '#fff3e0', text: '#E65100' },  // 橙色
        { bg: '#e8eaf6', border: '#e8eaf6', text: '#283593' },  // 靛蓝
        { bg: '#fce4ec', border: '#fce4ec', text: '#c2185b' },  // 粉红
        { bg: '#f1f8e9', border: '#f1f8e9', text: '#558b2f' }   // 浅绿
    ];

    // 定义需要显示的指标
    const desiredMetrics = {
        '中科院Top分区': 'sciUpTop',
        'JCR分区': 'sci',
        '5年影响因子': 'sciif5',
        '中科院基础版': 'sciBase',
        '中科院升级版': 'sciUp',
        '影响因子': 'sciif',
        'JCI指数': 'jci',
        'ESI学科分类': 'esi'
    };

    // 移除中文字符的函数
    function removeChineseCharacters(str) {
        return str.replace(/[\u4e00-\u9fa5]+/g, '').trim();
    }

    // 获取指标的className
    function getMetricClassName(metric) {
        return metric
            .toLowerCase()
            .replace(/[0-9]/g, '')
            .replace(/最新/g, '')
            .replace(/[：:]/g, '')
            .replace(/\s+/g, '-')
            .replace(/[年]/g, '')
            .trim();
    }

    // 添加缓存系统
    const CACHE_KEY = 'journal_cache';
    const CACHE_EXPIRY = 30 * 24 * 60 * 60 * 1000; // 30天过期时间(毫秒)
    const journalCache = new Map(loadCacheFromStorage());
    const requestQueue = [];
    let isProcessing = false;
    const RATE_LIMIT = 500; // 每次请求间隔500毫秒

    // 从localStorage加载缓存
    function loadCacheFromStorage() {
        try {
            const cached = localStorage.getItem(CACHE_KEY);
            if (!cached) return [];

            const parsedCache = JSON.parse(cached);
            const now = Date.now();

            // 过滤掉过期的数据
            const validEntries = Object.entries(parsedCache)
                .filter(([_, value]) => now - value.timestamp < CACHE_EXPIRY)
                .map(([key, value]) => [key, value.data]);

            // 清理过期数据
            if (validEntries.length !== Object.keys(parsedCache).length) {
                saveCacheToStorage(new Map(validEntries));
            }

            return validEntries;
        } catch (error) {
            console.error('加载缓存失败:', error);
            return [];
        }
    }

    // 保存缓存到localStorage
    function saveCacheToStorage(cache) {
        try {
            const cacheObject = {};
            cache.forEach((value, key) => {
                cacheObject[key] = {
                    data: value,
                    timestamp: Date.now()
                };
            });
            localStorage.setItem(CACHE_KEY, JSON.stringify(cacheObject));
        } catch (error) {
            console.error('保存缓存失败:', error);
        }
    }

    // 处理请求队列
    async function processQueue() {
        if (isProcessing || requestQueue.length === 0) return;

        isProcessing = true;
        while (requestQueue.length > 0) {
            const { journalName, resolve, reject } = requestQueue.shift();
            try {
                // 再次检查缓存(以防在队列中等待期间其他请求已经获取了数据)
                if (journalCache.has(journalName)) {
                    console.log('队列处理时从缓存获取:', journalName);
                    resolve(journalCache.get(journalName));
                    continue; // 跳过延迟直接处理下一个
                }

                const result = await executeRequest(journalName);
                resolve(result);
                // 只有实际发起API请求时才添加延迟
                await new Promise(resolve => setTimeout(resolve, RATE_LIMIT));
            } catch (error) {
                console.error('处理队列项失败:', error);
                reject(error);
            }
        }
        isProcessing = false;
    }

    // 执行实际的API请求
    async function executeRequest(journalName) {
        try {
            console.log('发起API请求:', journalName);
            const url = 'https://journal-api.yueyang.eu.org';
            const response = await fetch(url, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify({
                    journal: journalName
                })
            });

            if (!response.ok) {
                throw new Error(`API请求失败: ${response.status}`);
            }

            const data = await response.json();

            if (data.code === 200 && data.data && data.data.officialRank && data.data.officialRank.all) {
                const processedData = processJournalData(data.data.officialRank.all);
                if (processedData) {
                    console.log('成功获取数据:', journalName);
                    return processedData;
                }
            }
            console.log('无有效数据:', journalName);
            return null;
        } catch (error) {
            console.error('API请求失败:', journalName, error);
            return null;
        }
    }

    // 获取期刊数据
    async function getJournalData(journalName) {
        // 首先检查缓存
        if (journalCache.has(journalName)) {
            console.log('从缓存获取数据:', journalName);
            return journalCache.get(journalName);
        }

        return new Promise((resolve, reject) => {
            requestQueue.push({
                journalName,
                resolve: (data) => {
                    if (data) {
                        // 保存到缓存
                        journalCache.set(journalName, data);
                        saveCacheToStorage(journalCache);
                        console.log('数据已缓存:', journalName);
                    }
                    resolve(data);
                },
                reject
            });
            processQueue();
        });
    }

    // 处理期刊数据
    function processJournalData(metrics) {
        const result = {};
        for (const [displayName, key] of Object.entries(desiredMetrics)) {
            if (metrics[key] && metrics[key].toString().trim() !== '') {
                result[displayName] = metrics[key];
            }
        }
        return Object.keys(result).length > 0 ? result : null;
    }

    // 更新搜索页面的指标显示
    function updateSearchMetricsDisplay(container, journalData) {
        const style = document.createElement('style');
        style.textContent = `
            .journal-metrics {
                display: flex;
                flex-wrap: wrap;
                gap: 4px;
                margin-top: 4px;
                font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Arial, sans-serif;
            }
            .metric-item {
                padding: 2px 6px;
                border-radius: 12px;
                font-size: 12px;
                font-weight: 600;
                width: fit-content;
                display: flex;
                align-items: center;
            }
        `;

        Object.keys(desiredMetrics).forEach((metric, index) => {
            const className = getMetricClassName(metric);
            const colorScheme = colorSchemes[index % colorSchemes.length];
            style.textContent += `
                .metric-item.${className} {
                    background-color: ${colorScheme.bg};
                    border: none;
                    color: ${colorScheme.text};
                }
            `;
        });

        if (!document.head.querySelector('style[data-metrics-style]')) {
            style.setAttribute('data-metrics-style', 'true');
            document.head.appendChild(style);
        }

        const metricsContainer = document.createElement('div');
        metricsContainer.className = 'journal-metrics';

        Object.entries(journalData).forEach(([metric, value]) => {
            const metricElement = document.createElement('div');
            const className = getMetricClassName(metric);
            metricElement.className = `metric-item ${className}`;
            metricElement.textContent = `${metric}: ${value}`;
            metricsContainer.appendChild(metricElement);
        });

        container.appendChild(metricsContainer);
    }

    // 更新详情页面的指标显示
    function updateDetailMetricsDisplay(container, journalData) {
        const style = document.createElement('style');
        style.textContent = `
            .journal-metrics {
                display: flex;
                flex-wrap: wrap;
                gap: 8px;
                margin-top: 8px;
                font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Arial, sans-serif;
            }
            .metric-item {
                padding: 3px 6px;
                border-radius: 8px;
                font-size: 0.5em;
                width: fit-content;
                display: flex;
                align-items: center;
            }
        `;

        Object.keys(desiredMetrics).forEach((metric, index) => {
            const className = getMetricClassName(metric);
            const colorScheme = colorSchemes[index % colorSchemes.length];
            style.textContent += `
                .metric-item.${className} {
                    background-color: ${colorScheme.bg};
                    border: none;
                    color: ${colorScheme.text};
                }
            `;
        });

        if (!document.head.querySelector('style[data-metrics-style]')) {
            style.setAttribute('data-metrics-style', 'true');
            document.head.appendChild(style);
        }

        const metricsContainer = document.createElement('div');
        metricsContainer.className = 'journal-metrics';

        Object.entries(journalData).forEach(([metric, value]) => {
            const metricElement = document.createElement('div');
            const className = getMetricClassName(metric);
            metricElement.className = `metric-item ${className}`;
            metricElement.textContent = `${metric}: ${value}`;
            metricsContainer.appendChild(metricElement);
        });

        container.appendChild(metricsContainer);
    }

    // 修改搜索页面的显示逻辑
    async function addImpactFactorsToSearch() {
        console.log('执行搜索页面处理');
        const articles = document.querySelectorAll('.docsum-content');
        console.log('找到文章数量:', articles.length);

        for (const article of articles) {
            const citationSpan = article.querySelector('.docsum-journal-citation');
            if (!citationSpan ||
                citationSpan.querySelector('.journal-metrics') ||
                citationSpan.dataset.processing === 'true') continue;

            const citationText = citationSpan.textContent;
            const journalRegex = /^([^.]+)\./;
            const match = citationText.match(journalRegex);

            if (match && match[1]) {
                const journalName = removeChineseCharacters(match[1]);
                console.log('处理期刊:', journalName);
                citationSpan.dataset.processing = 'true';

                try {
                    if (journalCache.has(journalName)) {
                        const journalData = journalCache.get(journalName);
                        if (journalData) {
                            updateSearchMetricsDisplay(citationSpan, journalData);
                        }
                        delete citationSpan.dataset.processing;
                        continue;
                    }

                    const journalData = await getJournalData(journalName);
                    if (journalData) {
                        updateSearchMetricsDisplay(citationSpan, journalData);
                    }
                } catch (error) {
                    console.error('处理期刊数据失败:', error);
                } finally {
                    delete citationSpan.dataset.processing;
                }
            }
        }
    }

    // 修改详情页面的显示逻辑
    async function addImpactFactorsToDetail() {
        console.log('执行详情页面处理');
        const journalElement = document.getElementById('full-view-journal-trigger');
        const existingMetrics = document.querySelector('.journal-metrics');

        if (!journalElement ||
            existingMetrics ||
            journalElement.dataset.processing === 'true') return;

        const journalName = removeChineseCharacters(journalElement.textContent);
        console.log('处理期刊:', journalName);
        journalElement.dataset.processing = 'true';

        try {
            const journalData = await getJournalData(journalName);
            if (journalData) {
                const titleElement = document.querySelector('.heading-title');
                if (titleElement) {
                    updateDetailMetricsDisplay(titleElement, journalData);
                }
            }
        } catch (error) {
            console.error('处理期刊数据失败:', error);
        } finally {
            delete journalElement.dataset.processing;
        }
    }

    // 添加PMC页面的期刊信息处理
    async function addImpactFactorsToPMC() {
        console.log('执行PMC页面处理');
        const journalElement = document.querySelector('main article section:first-child div button');
        const titleContainer = document.querySelector('main article section:nth-child(2) div hgroup');

        if (!journalElement || !titleContainer || titleContainer.querySelector('.journal-metrics')) {
            return;
        }

        const journalName = removeChineseCharacters(journalElement.textContent);
        console.log('处理期刊:', journalName);

        try {
            const journalData = await getJournalData(journalName);
            if (journalData) {
                const style = document.createElement('style');
                style.textContent = `
                    .journal-metrics {
                        display: flex;
                        flex-wrap: wrap;
                        gap: 8px;
                        margin-top: 8px;
                        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Arial, sans-serif;
                    }
                    .metric-item {
                        padding: 3px 8px;
                        border-radius: 8px;
                        font-size: 0.9em;
                        width: fit-content;
                        display: flex;
                        align-items: center;
                    }
                `;

                Object.keys(desiredMetrics).forEach((metric, index) => {
                    const className = getMetricClassName(metric);
                    const colorScheme = colorSchemes[index % colorSchemes.length];
                    style.textContent += `
                        .metric-item.${className} {
                            background-color: ${colorScheme.bg};
                            border: none;
                            color: ${colorScheme.text};
                        }
                    `;
                });

                if (!document.head.querySelector('style[data-metrics-style]')) {
                    style.setAttribute('data-metrics-style', 'true');
                    document.head.appendChild(style);
                }

                const metricsContainer = document.createElement('div');
                metricsContainer.className = 'journal-metrics';

                Object.entries(journalData).forEach(([metric, value]) => {
                    const metricElement = document.createElement('div');
                    const className = getMetricClassName(metric);
                    metricElement.className = `metric-item ${className}`;
                    metricElement.textContent = `${metric}: ${value}`;
                    metricsContainer.appendChild(metricElement);
                });

                titleContainer.appendChild(metricsContainer);
            }
        } catch (error) {
            console.error('处理PMC期刊数据失败:', error);
        }
    }

    // 主函数
    async function init() {
        console.log('init 函数开始执行');
        const isPMC = window.location.hostname.includes('pmc.ncbi.nlm.nih.gov');
        const isPubMedDetail = window.location.pathname.match(/\/\d+\//);
        console.log('页面类型:', { isPMC, isPubMedDetail, location: window.location.href });

        if (isPMC) {
            await addImpactFactorsToPMC();
        } else if (isPubMedDetail) {
            await addImpactFactorsToDetail();
        } else {
            await addImpactFactorsToSearch();
        }
    }

    // 在页面加载完成后执行
    document.addEventListener('DOMContentLoaded', () => {
        console.log('DOMContentLoaded 触发');
        init();
    });

    window.addEventListener('load', () => {
        console.log('Window load 触发');
        init();
    });

    // 创建观察器处理动态加载的内容
    const observer = new MutationObserver((mutations) => {
        console.log('检测到页面变化');
        setTimeout(init, 500);
    });

    // 确保 observer 正确设置
    setTimeout(() => {
        const target = document.querySelector('.search-results-chunk') ||
                      document.querySelector('main') ||
                      document.body;

        if (target) {
            console.log('开始观察DOM变化');
            observer.observe(target, {
                childList: true,
                subtree: true,
                attributes: true
            });
        } else {
            console.log('未找到观察目标');
        }
    }, 1000);
})();
